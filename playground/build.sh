#!/bin/sh
# Builds the playground as static files in a directory, target/playground by default.
#
#     playground/build.sh [out]
set -eu

root=$(cd "$(dirname "$0")/.." && pwd)
target=${CARGO_TARGET_DIR:-$root/target}
out=${1:-$target/playground}

cargo build --manifest-path "$root/Cargo.toml" -p moiety-wasm --target-dir "$target" \
  --target wasm32-unknown-unknown --profile wasm
mkdir -p "$out"
cp "$target/wasm32-unknown-unknown/wasm/moiety_wasm.wasm" "$out/moiety.wasm"
cp "$root/playground/index.html" "$root/playground/playground.css" \
  "$root/playground/playground.js" "$root/playground/worker.js" "$out/"
cp "$root/crates/moiety-draw/src/html/viewer.css" \
  "$root/crates/moiety-draw/src/html/viewer.js" "$out/"
