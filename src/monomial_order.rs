//! Monomial orderings: lex, grlex, grevlex.

use serde::{Deserialize, Serialize};
use std::cmp::Ordering;

/// Supported monomial orderings.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MonomialOrder {
    /// Lexicographic order: compares variables left-to-right.
    Lex,
    /// Graded lexicographic order: total degree first, then lex.
    GrLex,
    /// Graded reverse lexicographic order: total degree first, then reverse lex.
    GrRevLex,
}

impl MonomialOrder {
    /// Compare two exponent vectors under this ordering.
    pub fn compare(&self, a: &[u32], b: &[u32]) -> Ordering {
        debug_assert_eq!(a.len(), b.len(), "exponent vectors must have same length");
        match self {
            MonomialOrder::Lex => {
                for i in 0..a.len() {
                    match a[i].cmp(&b[i]) {
                        Ordering::Equal => continue,
                        o => return o,
                    }
                }
                Ordering::Equal
            }
            MonomialOrder::GrLex => {
                let deg_a: u32 = a.iter().sum();
                let deg_b: u32 = b.iter().sum();
                match deg_a.cmp(&deg_b) {
                    Ordering::Equal => {}
                    o => return o,
                }
                // tiebreak with lex
                for i in 0..a.len() {
                    match a[i].cmp(&b[i]) {
                        Ordering::Equal => continue,
                        o => return o,
                    }
                }
                Ordering::Equal
            }
            MonomialOrder::GrRevLex => {
                let deg_a: u32 = a.iter().sum();
                let deg_b: u32 = b.iter().sum();
                match deg_a.cmp(&deg_b) {
                    Ordering::Equal => {}
                    o => return o,
                }
                // tiebreak: prefer smaller exponent in the LAST variable
                for i in (0..a.len()).rev() {
                    match b[i].cmp(&a[i]) {
                        Ordering::Equal => continue,
                        o => return o,
                    }
                }
                Ordering::Equal
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn exponents(e: &[u32]) -> std::borrow::Cow<'_, [u32]> {
        std::borrow::Cow::Borrowed(e)
    }

    #[test]
    fn test_lex_basic() {
        // x > y: [1,0] > [0,1] under lex
        assert_eq!(MonomialOrder::Lex.compare(&[1, 0], &[0, 1]), Ordering::Greater);
        // x^2 > x*y
        assert_eq!(MonomialOrder::Lex.compare(&[2, 0], &[1, 1]), Ordering::Greater);
        // equal
        assert_eq!(MonomialOrder::Lex.compare(&[1, 2], &[1, 2]), Ordering::Equal);
    }

    #[test]
    fn test_grlex_basic() {
        // same degree, tiebreak lex: x^2 > x*y
        assert_eq!(MonomialOrder::GrLex.compare(&[2, 0], &[1, 1]), Ordering::Greater);
        // higher degree wins: x*y (deg 2) > x (deg 1)
        assert_eq!(MonomialOrder::GrLex.compare(&[1, 1], &[1, 0]), Ordering::Greater);
        assert_eq!(MonomialOrder::GrLex.compare(&[1, 0], &[1, 1]), Ordering::Less);
    }

    #[test]
    fn test_grevlex_basic() {
        // same degree (2), grevlex: prefer smaller last variable exponent
        // x^2 [2,0] vs x*y [1,1]: last exponent 0 < 1, so x^2 > x*y in grevlex
        assert_eq!(MonomialOrder::GrRevLex.compare(&[2, 0], &[1, 1]), Ordering::Greater);
        // x*y [1,1] vs y^2 [0,2]: last exponent 1 < 2, so x*y > y^2
        assert_eq!(MonomialOrder::GrRevLex.compare(&[1, 1], &[0, 2]), Ordering::Greater);
        // higher degree always wins
        assert_eq!(MonomialOrder::GrRevLex.compare(&[1, 1], &[1, 0]), Ordering::Greater);
    }

    #[test]
    fn test_lex_three_vars() {
        // x > y > z: [1,0,0] > [0,1,0] > [0,0,1]
        assert_eq!(MonomialOrder::Lex.compare(&[1, 0, 0], &[0, 1, 0]), Ordering::Greater);
        assert_eq!(MonomialOrder::Lex.compare(&[0, 1, 0], &[0, 0, 1]), Ordering::Greater);
    }

    #[test]
    fn test_zero_exponents() {
        assert_eq!(MonomialOrder::Lex.compare(&[0, 0], &[0, 0]), Ordering::Equal);
        assert_eq!(MonomialOrder::GrLex.compare(&[0, 0], &[0, 0]), Ordering::Equal);
        assert_eq!(MonomialOrder::GrRevLex.compare(&[0, 0], &[0, 0]), Ordering::Equal);
    }
}
