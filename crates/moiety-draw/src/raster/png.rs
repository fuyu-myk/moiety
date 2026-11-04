//! PNG encoding of RGBA pixels, compressed with fixed-code DEFLATE.
//!
//! Matches are found greedily through a table of the last position of every
//! three-byte prefix, which suits line drawings: long runs of one colour
//! shrink to a few codes.

/// Shortest and longest match DEFLATE can code.
const MIN_MATCH: usize = 3;
const MAX_MATCH: usize = 258;

/// Furthest back a match may start.
const WINDOW: usize = 32_768;

/// Bits of the three-byte prefix hash.
const HASH_BITS: u32 = 15;

/// Base lengths of the length codes 257 to 285, and their extra bits.
const LENGTH_BASE: [u16; 29] = [
    3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131,
    163, 195, 227, 258,
];
const LENGTH_EXTRA: [u8; 29] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
];

/// Base distances of the distance codes 0 to 29, and their extra bits.
const DISTANCE_BASE: [u16; 30] = [
    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537,
    2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
];
const DISTANCE_EXTRA: [u8; 30] = [
    0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13,
    13,
];

/// A PNG file of the RGBA `pixels`, `width` pixels to a row.
pub(super) fn encode(width: u32, height: u32, pixels: &[u8]) -> Vec<u8> {
    let row = width as usize * 4;
    let mut raw = Vec::with_capacity((row + 1) * height as usize);
    for line in pixels.chunks(row.max(1)).take(height as usize) {
        raw.push(0);
        raw.extend_from_slice(line);
    }
    let mut header = Vec::with_capacity(13);
    header.extend_from_slice(&width.to_be_bytes());
    header.extend_from_slice(&height.to_be_bytes());
    header.extend_from_slice(&[8, 6, 0, 0, 0]);
    let mut png = b"\x89PNG\r\n\x1a\n".to_vec();
    chunk(&mut png, b"IHDR", &header);
    chunk(&mut png, b"IDAT", &zlib(&raw));
    chunk(&mut png, b"IEND", &[]);
    png
}

fn chunk(png: &mut Vec<u8>, kind: &[u8; 4], data: &[u8]) {
    png.extend_from_slice(&(data.len() as u32).to_be_bytes());
    let start = png.len();
    png.extend_from_slice(kind);
    png.extend_from_slice(data);
    let crc = crc32(&png[start..]);
    png.extend_from_slice(&crc.to_be_bytes());
}

/// `data` as a zlib stream.
fn zlib(data: &[u8]) -> Vec<u8> {
    let mut out = vec![0x78, 0x01];
    out.extend(deflate(data));
    out.extend_from_slice(&adler32(data).to_be_bytes());
    out
}

/// `data` as one final DEFLATE block with the fixed codes.
fn deflate(data: &[u8]) -> Vec<u8> {
    let mut bits = Bits::default();
    bits.push(0b011, 3);
    let mut last = vec![usize::MAX; 1 << HASH_BITS];
    let mut at = 0;
    while at < data.len() {
        let found = (at + MIN_MATCH <= data.len()).then(|| {
            let key = hash(&data[at..at + MIN_MATCH]);
            let candidate = last[key];
            last[key] = at;
            candidate
        });
        let length = found
            .filter(|&from| from != usize::MAX && at - from <= WINDOW)
            .map_or(0, |from| {
                let limit = (data.len() - at).min(MAX_MATCH);
                (0..limit)
                    .take_while(|&k| data[from + k] == data[at + k])
                    .count()
                    .min(limit)
            });
        if length >= MIN_MATCH {
            let distance = at - found.unwrap_or(at);
            bits.length(length);
            bits.distance(distance);
            for skipped in at + 1..(at + length).min(data.len().saturating_sub(MIN_MATCH - 1)) {
                last[hash(&data[skipped..skipped + MIN_MATCH])] = skipped;
            }
            at += length;
        } else {
            bits.symbol(u16::from(data[at]));
            at += 1;
        }
    }
    bits.symbol(256);
    bits.finish()
}

fn hash(bytes: &[u8]) -> usize {
    let word = u32::from(bytes[0]) << 16 | u32::from(bytes[1]) << 8 | u32::from(bytes[2]);
    (word.wrapping_mul(2_654_435_761) >> (32 - HASH_BITS)) as usize
}

/// A DEFLATE bit stream, filled from the least significant bit of each byte.
#[derive(Default)]
struct Bits {
    bytes: Vec<u8>,
    buffer: u64,
    count: u32,
}

impl Bits {
    fn push(&mut self, value: u32, count: u32) {
        self.buffer |= u64::from(value) << self.count;
        self.count += count;
        while self.count >= 8 {
            self.bytes.push(self.buffer as u8);
            self.buffer >>= 8;
            self.count -= 8;
        }
    }

    /// A Huffman code, which DEFLATE stores from its most significant bit.
    fn code(&mut self, code: u32, count: u32) {
        self.push(code.reverse_bits() >> (32 - count), count);
    }

    /// A literal byte, or the end-of-block or a length symbol, in the fixed code.
    fn symbol(&mut self, symbol: u16) {
        let symbol = u32::from(symbol);
        match symbol {
            0..=143 => self.code(0x30 + symbol, 8),
            144..=255 => self.code(0x190 + symbol - 144, 9),
            256..=279 => self.code(symbol - 256, 7),
            _ => self.code(0xc0 + symbol - 280, 8),
        }
    }

    fn length(&mut self, length: usize) {
        let index = LENGTH_BASE.partition_point(|&base| usize::from(base) <= length) - 1;
        self.symbol(257 + index as u16);
        let extra = u32::from(LENGTH_EXTRA[index]);
        self.push((length - usize::from(LENGTH_BASE[index])) as u32, extra);
    }

    fn distance(&mut self, distance: usize) {
        let index = DISTANCE_BASE.partition_point(|&base| usize::from(base) <= distance) - 1;
        self.code(index as u32, 5);
        let extra = u32::from(DISTANCE_EXTRA[index]);
        self.push((distance - usize::from(DISTANCE_BASE[index])) as u32, extra);
    }

    fn finish(mut self) -> Vec<u8> {
        if self.count > 0 {
            self.bytes.push(self.buffer as u8);
        }
        self.bytes
    }
}

fn crc32(bytes: &[u8]) -> u32 {
    !bytes.iter().fold(u32::MAX, |crc, &byte| {
        (0..8).fold(crc ^ u32::from(byte), |crc, _| {
            (crc >> 1) ^ (0xedb8_8320 & (crc & 1).wrapping_neg())
        })
    })
}

fn adler32(bytes: &[u8]) -> u32 {
    const MOD: u32 = 65_521;
    let (a, b) = bytes
        .chunks(5552)
        .fold((1u32, 0u32), |(mut a, mut b), block| {
            for &byte in block {
                a += u32::from(byte);
                b += a;
            }
            (a % MOD, b % MOD)
        });
    b << 16 | a
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reads a stream of fixed-code DEFLATE blocks.
    fn inflate(stream: &[u8]) -> Vec<u8> {
        let mut at = 0;
        let mut bit = |count: u32| {
            (0..count).fold(0u32, |value, k| {
                let b = (stream[at / 8] >> (at % 8)) & 1;
                at += 1;
                value | u32::from(b) << k
            })
        };
        let mut out = Vec::new();
        loop {
            let last = bit(1);
            assert_eq!(bit(2), 1, "fixed codes only");
            loop {
                let mut code = 0;
                let mut length = 0;
                let symbol = loop {
                    code = code << 1 | bit(1);
                    length += 1;
                    match (length, code) {
                        (7, 0..=0x17) => break code + 256,
                        (8, 0x30..=0xbf) => break code - 0x30,
                        (8, 0xc0..=0xc7) => break code - 0xc0 + 280,
                        (9, 0x190..=0x1ff) => break code - 0x190 + 144,
                        _ => {}
                    }
                };
                match symbol {
                    0..=255 => out.push(symbol as u8),
                    256 => break,
                    _ => {
                        let i = (symbol - 257) as usize;
                        let length =
                            usize::from(LENGTH_BASE[i]) + bit(u32::from(LENGTH_EXTRA[i])) as usize;
                        let d = (0..5).fold(0, |v, _| v << 1 | bit(1)) as usize;
                        let distance = usize::from(DISTANCE_BASE[d])
                            + bit(u32::from(DISTANCE_EXTRA[d])) as usize;
                        for _ in 0..length {
                            out.push(out[out.len() - distance]);
                        }
                    }
                }
            }
            if last == 1 {
                return out;
            }
        }
    }

    #[test]
    fn checksums_match_their_reference_values() {
        assert_eq!(crc32(b"123456789"), 0xcbf4_3926);
        assert_eq!(crc32(b"IEND"), 0xae42_6082);
        assert_eq!(adler32(b"Wikipedia"), 0x11e6_0398);
        assert_eq!(adler32(&[0xff; 100_000]), {
            let (a, b) = (0..100_000u64).fold((1u64, 0u64), |(a, b), _| {
                ((a + 255) % 65_521, (b + a + 255) % 65_521)
            });
            (b << 16 | a) as u32
        });
    }

    #[test]
    fn deflates_and_inflates_back() {
        let text = b"a line drawing is mostly background, background, background".repeat(20);
        let mut noise: Vec<u8> = (0..70_000u32)
            .map(|i| (i.wrapping_mul(2_654_435_761) >> 24) as u8)
            .collect();
        noise.extend(vec![0; 1000]);
        for data in [
            &b""[..],
            b"x",
            b"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
            &text,
            &noise,
        ] {
            let packed = deflate(data);
            assert_eq!(inflate(&packed), data);
        }
        assert!(deflate(&text).len() < text.len() / 5);
        assert!(deflate(&[7; 100_000]).len() < 1000);
    }

    #[test]
    fn writes_a_well_formed_png() {
        let pixels = [255, 0, 0, 255, 0, 0, 0, 0].repeat(3);
        let png = encode(2, 3, &pixels);
        assert!(png.starts_with(b"\x89PNG\r\n\x1a\n\0\0\0\x0dIHDR\0\0\0\x02\0\0\0\x03\x08\x06"));
        assert!(png.ends_with(b"\0\0\0\0IEND\xae\x42\x60\x82"));
        let length = u32::from_be_bytes(png[33..37].try_into().unwrap()) as usize;
        assert_eq!(&png[37..41], b"IDAT");
        let stream = &png[41..41 + length];
        assert_eq!(&stream[..2], &[0x78, 0x01]);
        let raw = inflate(&stream[2..length - 4]);
        assert_eq!(raw, [&[0][..], &pixels[..8]].concat().repeat(3));
        assert_eq!(stream[length - 4..], adler32(&raw).to_be_bytes());
        let crc = crc32(&png[37..41 + length]);
        assert_eq!(png[41 + length..45 + length], crc.to_be_bytes());
    }
}
