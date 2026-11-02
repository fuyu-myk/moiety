use std::fmt;

/// Broad bonding behaviour of an element, used by the valence model.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Category {
    /// Nonmetals, metalloids and noble gases.
    Nonmetal,
    /// s-block and p-block metals, plus zinc, cadmium and mercury. These may bond to carbon.
    MainGroupMetal,
    /// Groups 3 to 11.
    TransitionMetal,
    /// Lanthanides and actinides.
    FBlock,
    /// Atomic number 104 and above.
    Superheavy,
}

struct Data {
    symbol: &'static str,
    valences: &'static [u8],
    category: Category,
}

const fn d(symbol: &'static str, valences: &'static [u8], category: Category) -> Data {
    Data {
        symbol,
        valences,
        category,
    }
}

use Category::{
    FBlock as F, MainGroupMetal as M, Nonmetal as N, Superheavy as X, TransitionMetal as T,
};

static TABLE: [Data; 118] = [
    d("H", &[1], N),
    d("He", &[0], N),
    d("Li", &[1], M),
    d("Be", &[2], M),
    d("B", &[3], N),
    d("C", &[4], N),
    d("N", &[3, 5], N),
    d("O", &[2], N),
    d("F", &[1], N),
    d("Ne", &[0], N),
    d("Na", &[1], M),
    d("Mg", &[2], M),
    d("Al", &[3], M),
    d("Si", &[4], N),
    d("P", &[3, 5], N),
    d("S", &[2, 4, 6], N),
    d("Cl", &[1], N),
    d("Ar", &[0], N),
    d("K", &[1], M),
    d("Ca", &[2], M),
    d("Sc", &[3], T),
    d("Ti", &[2, 3, 4], T),
    d("V", &[2, 3, 4, 5], T),
    d("Cr", &[2, 3, 6], T),
    d("Mn", &[2, 3, 4, 7], T),
    d("Fe", &[2, 3], T),
    d("Co", &[2, 3], T),
    d("Ni", &[2], T),
    d("Cu", &[1, 2], T),
    d("Zn", &[2], M),
    d("Ga", &[3], M),
    d("Ge", &[4], N),
    d("As", &[3, 5], N),
    d("Se", &[2, 4, 6], N),
    d("Br", &[1], N),
    d("Kr", &[0], N),
    d("Rb", &[1], M),
    d("Sr", &[2], M),
    d("Y", &[3], T),
    d("Zr", &[4], T),
    d("Nb", &[3, 5], T),
    d("Mo", &[4, 6], T),
    d("Tc", &[4, 7], T),
    d("Ru", &[2, 3, 4], T),
    d("Rh", &[3], T),
    d("Pd", &[2, 4], T),
    d("Ag", &[1], T),
    d("Cd", &[2], M),
    d("In", &[3], M),
    d("Sn", &[2, 4], M),
    d("Sb", &[3, 5], N),
    d("Te", &[2, 4, 6], N),
    d("I", &[1, 3, 5], N),
    d("Xe", &[0, 2, 4, 6], N),
    d("Cs", &[1], M),
    d("Ba", &[2], M),
    d("La", &[3], F),
    d("Ce", &[3, 4], F),
    d("Pr", &[3], F),
    d("Nd", &[3], F),
    d("Pm", &[3], F),
    d("Sm", &[2, 3], F),
    d("Eu", &[2, 3], F),
    d("Gd", &[3], F),
    d("Tb", &[3], F),
    d("Dy", &[3], F),
    d("Ho", &[3], F),
    d("Er", &[3], F),
    d("Tm", &[3], F),
    d("Yb", &[2, 3], F),
    d("Lu", &[3], F),
    d("Hf", &[4], T),
    d("Ta", &[5], T),
    d("W", &[4, 6], T),
    d("Re", &[4, 7], T),
    d("Os", &[4, 8], T),
    d("Ir", &[3, 4], T),
    d("Pt", &[2, 4], T),
    d("Au", &[1, 3], T),
    d("Hg", &[1, 2], M),
    d("Tl", &[1, 3], M),
    d("Pb", &[2, 4], M),
    d("Bi", &[3, 5], M),
    d("Po", &[2, 4], M),
    d("At", &[1], N),
    d("Rn", &[0], N),
    d("Fr", &[1], M),
    d("Ra", &[2], M),
    d("Ac", &[3], F),
    d("Th", &[4], F),
    d("Pa", &[5], F),
    d("U", &[3, 4, 6], F),
    d("Np", &[5], F),
    d("Pu", &[4], F),
    d("Am", &[3], F),
    d("Cm", &[3], F),
    d("Bk", &[3, 4], F),
    d("Cf", &[3], F),
    d("Es", &[3], F),
    d("Fm", &[3], F),
    d("Md", &[3], F),
    d("No", &[2, 3], F),
    d("Lr", &[3], F),
    d("Rf", &[4], X),
    d("Db", &[5], X),
    d("Sg", &[6], X),
    d("Bh", &[7], X),
    d("Hs", &[8], X),
    d("Mt", &[3], X),
    d("Ds", &[2], X),
    d("Rg", &[1], X),
    d("Cn", &[2], X),
    d("Nh", &[1, 3], X),
    d("Fl", &[2, 4], X),
    d("Mc", &[1, 3], X),
    d("Lv", &[2, 4], X),
    d("Ts", &[1], X),
    d("Og", &[0], X),
];

/// A chemical element, identified by its atomic number.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Element(u8);

impl Element {
    /// Hydrogen.
    pub const H: Element = Element(1);
    /// Boron.
    pub const B: Element = Element(5);
    /// Carbon.
    pub const C: Element = Element(6);
    /// Nitrogen.
    pub const N: Element = Element(7);
    /// Oxygen.
    pub const O: Element = Element(8);
    /// Fluorine.
    pub const F: Element = Element(9);
    /// Silicon.
    pub const SI: Element = Element(14);
    /// Phosphorus.
    pub const P: Element = Element(15);
    /// Sulfur.
    pub const S: Element = Element(16);
    /// Chlorine.
    pub const CL: Element = Element(17);
    /// Selenium.
    pub const SE: Element = Element(34);
    /// Bromine.
    pub const BR: Element = Element(35);
    /// Iodine.
    pub const I: Element = Element(53);

    /// The element with atomic number `z`, if `z` is between 1 and 118.
    #[must_use]
    pub fn from_atomic_number(z: u8) -> Option<Element> {
        (1..=118).contains(&z).then_some(Element(z))
    }

    /// The element whose symbol is exactly `symbol`, with standard capitalization.
    #[must_use]
    pub fn from_symbol(symbol: &str) -> Option<Element> {
        TABLE
            .iter()
            .position(|data| data.symbol == symbol)
            .and_then(|i| u8::try_from(i + 1).ok())
            .map(Element)
    }

    /// Atomic number.
    #[must_use]
    pub const fn atomic_number(self) -> u8 {
        self.0
    }

    /// Symbol with standard capitalization, such as `Cl`.
    #[must_use]
    pub fn symbol(self) -> &'static str {
        self.data().symbol
    }

    /// Bonding category.
    #[must_use]
    pub fn category(self) -> Category {
        self.data().category
    }

    /// Allowed valences of the neutral atom, lowest first.
    #[must_use]
    pub fn valences(self) -> &'static [u8] {
        self.data().valences
    }

    /// Allowed valences of the atom carrying formal charge `charge`, lowest first.
    ///
    /// Nonmetals take the lowest valence of the isoelectronic element `Z - charge`.
    /// Metals lose one valence per unit of positive charge and take no negative charge.
    #[must_use]
    pub fn charged_valences(self, charge: i8) -> Vec<u8> {
        if charge == 0 {
            return self.valences().to_vec();
        }
        if self.is_metal() {
            let Ok(loss) = u8::try_from(charge) else {
                return Vec::new();
            };
            return self
                .valences()
                .iter()
                .filter_map(|v| v.checked_sub(loss))
                .collect();
        }
        let shifted = i16::from(self.0) - i16::from(charge);
        match u8::try_from(shifted)
            .ok()
            .and_then(Element::from_atomic_number)
        {
            Some(iso) => vec![iso.valences()[0]],
            None => vec![0],
        }
    }

    /// Whether the element is one of H B C N O F Si P S Cl Br I.
    #[must_use]
    pub fn is_organic(self) -> bool {
        matches!(self.0, 1 | 5 | 6 | 7 | 8 | 9 | 14 | 15 | 16 | 17 | 35 | 53)
    }

    /// Whether the element is a metal of any kind.
    #[must_use]
    pub fn is_metal(self) -> bool {
        self.category() != Category::Nonmetal
    }

    /// Whether the element can take part in double and triple bonds.
    #[must_use]
    pub fn allows_multiple_bonds(self) -> bool {
        matches!(
            self.0,
            5 | 6 | 7 | 8 | 14 | 15 | 16 | 17 | 32 | 33 | 34 | 35 | 51 | 52 | 53 | 54
        )
    }

    /// Whether the element may bond to carbon and hydrogen.
    ///
    /// Transition metals, f-block elements and superheavy elements bond only to heteroatoms.
    #[must_use]
    pub fn bonds_to_carbon(self) -> bool {
        matches!(
            self.category(),
            Category::Nonmetal | Category::MainGroupMetal
        )
    }

    fn data(self) -> &'static Data {
        &TABLE[usize::from(self.0) - 1]
    }
}

impl fmt::Debug for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.symbol())
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.symbol())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symbols_round_trip() {
        for z in 1..=118 {
            let element = Element::from_atomic_number(z).unwrap();
            assert_eq!(Element::from_symbol(element.symbol()), Some(element));
        }
        assert_eq!(Element::from_symbol("CL"), None);
        assert_eq!(Element::from_atomic_number(0), None);
    }

    #[test]
    fn charged_valences_follow_isoelectronic_neighbours() {
        assert_eq!(Element::N.charged_valences(1), vec![4]);
        assert_eq!(Element::O.charged_valences(-1), vec![1]);
        assert_eq!(Element::C.charged_valences(-1), vec![3]);
        assert_eq!(Element::B.charged_valences(-1), vec![4]);
        assert_eq!(Element::CL.charged_valences(-1), vec![0]);
        assert_eq!(Element::H.charged_valences(1), vec![0]);
    }

    #[test]
    fn metals_lose_valence_with_charge() {
        let iron = Element::from_symbol("Fe").unwrap();
        assert_eq!(iron.charged_valences(3), vec![0]);
        assert_eq!(iron.charged_valences(2), vec![0, 1]);
        assert!(iron.charged_valences(-1).is_empty());
    }

    #[test]
    fn transition_metals_avoid_carbon() {
        assert!(!Element::from_symbol("Co").unwrap().bonds_to_carbon());
        assert!(Element::from_symbol("Mg").unwrap().bonds_to_carbon());
        assert!(Element::from_symbol("Zn").unwrap().bonds_to_carbon());
        assert!(!Element::from_symbol("Ts").unwrap().bonds_to_carbon());
    }
}
