//! Polynomial ideals and ideal membership testing.

use serde::{Deserialize, Serialize};
use crate::polynomial::{Polynomial, PolynomialRing, Term};
use crate::monomial_order::MonomialOrder;
use crate::groebner::reduce_by_set;

/// A polynomial ideal.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ideal {
    /// Generators of the ideal.
    pub generators: Vec<Polynomial>,
    /// The polynomial ring this ideal lives in.
    pub ring: PolynomialRing,
}

impl Ideal {
    /// Create a new ideal from generators.
    pub fn new(generators: Vec<Polynomial>, ring: PolynomialRing) -> Self {
        let gens: Vec<Polynomial> = generators
            .into_iter()
            .filter(|p| !p.is_zero())
            .collect();
        Self { generators: gens, ring }
    }

    /// The zero ideal {0}.
    pub fn zero(ring: PolynomialRing) -> Self {
        Self {
            generators: vec![],
            ring,
        }
    }

    /// Check if the ideal is the zero ideal.
    pub fn is_zero(&self) -> bool {
        self.generators.is_empty()
    }

    /// Number of generators.
    pub fn num_generators(&self) -> usize {
        self.generators.len()
    }

    /// Check ideal membership using Groebner basis reduction.
    /// A polynomial f is in the ideal I iff its normal form w.r.t. a Groebner basis of I is 0.
    pub fn contains(&self, f: &Polynomial) -> bool {
        if self.is_zero() {
            return f.is_zero();
        }
        use crate::groebner::groebner_basis;
        let gb = groebner_basis(&self.generators, self.ring.order);
        let (_, remainder) = reduce_by_set(f, &gb);
        remainder.is_zero()
    }

    /// Sum of two ideals: I + J = <generators of I ∪ generators of J>.
    pub fn sum(&self, other: &Ideal) -> Ideal {
        debug_assert_eq!(self.ring.num_vars, other.ring.num_vars);
        let mut gens = self.generators.clone();
        gens.extend(other.generators.iter().cloned());
        Ideal::new(gens, self.ring.clone())
    }

    /// Product of two ideals: I * J.
    pub fn product(&self, other: &Ideal) -> Ideal {
        debug_assert_eq!(self.ring.num_vars, other.ring.num_vars);
        let mut gens = Vec::new();
        for a in &self.generators {
            for b in &other.generators {
                gens.push(a.mul(b));
            }
        }
        Ideal::new(gens, self.ring.clone())
    }

    /// Intersection of two ideals using Groebner basis elimination.
    /// I ∩ J = (t*I, (1-t)*J) ∩ k[x1,...,xn] eliminating t.
    pub fn intersection(&self, other: &Ideal) -> Ideal {
        debug_assert_eq!(self.ring.num_vars, other.ring.num_vars);
        let n = self.ring.num_vars;
        // Use auxiliary variable t at index n
        let mut order = self.ring.order;
        let aux_ring = PolynomialRing::new(n + 1, order);

        // t*I: multiply each generator by t
        let mut augmented = Vec::new();
        let t = Polynomial::variable(n, n + 1, order);
        for g in &self.generators {
            // extend g to n+1 vars
            let mut ext_g = g.clone();
            ext_g.num_vars = n + 1;
            for term in &mut ext_g.terms {
                term.exponents.push(0);
            }
            augmented.push(t.mul(&ext_g));
        }

        // (1-t)*J
        let one = Polynomial::constant(1.0, n + 1, order);
        let one_minus_t = one.sub(&t);
        for g in &other.generators {
            let mut ext_g = g.clone();
            ext_g.num_vars = n + 1;
            for term in &mut ext_g.terms {
                term.exponents.push(0);
            }
            augmented.push(one_minus_t.mul(&ext_g));
        }

        // Compute Groebner basis and eliminate t
        use crate::groebner::groebner_basis;
        let gb = groebner_basis(&augmented, order);

        // Select elements not involving t (exponent at index n is 0)
        let result_gens: Vec<Polynomial> = gb
            .into_iter()
            .filter(|p| {
                p.terms.iter().all(|t| t.exponents.get(n).copied().unwrap_or(0) == 0)
            })
            .map(|mut p| {
                // drop the auxiliary variable
                p.num_vars = n;
                for term in &mut p.terms {
                    term.exponents.pop();
                }
                p
            })
            .collect();

        Ideal::new(result_gens, self.ring.clone())
    }

    /// Check if this ideal contains another ideal (I ⊇ J).
    pub fn contains_ideal(&self, other: &Ideal) -> bool {
        other.generators.iter().all(|g| self.contains(g))
    }

    /// Check if two ideals are equal.
    pub fn equals(&self, other: &Ideal) -> bool {
        self.contains_ideal(other) && other.contains_ideal(self)
    }

    /// The radical ideal test (simplified): check if f^m in I for some m.
    /// This is a simplified check — just tests f in I for basic radical membership.
    pub fn radical_contains(&self, f: &Polynomial, max_power: u32) -> bool {
        let mut power = f.clone();
        for _ in 0..max_power {
            if self.contains(&power) {
                return true;
            }
            power = power.mul(f);
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::monomial_order::MonomialOrder;

    fn ring2() -> PolynomialRing {
        PolynomialRing::new(2, MonomialOrder::Lex)
    }

    #[test]
    fn test_ideal_creation() {
        let ring = ring2();
        let f = ring.from_raw(vec![(1.0, vec![2, 0]), (-1.0, vec![0, 0])]); // x^2 - 1
        let ideal = Ideal::new(vec![f], ring);
        assert_eq!(ideal.num_generators(), 1);
    }

    #[test]
    fn test_zero_ideal() {
        let ring = ring2();
        let ideal = Ideal::zero(ring);
        assert!(ideal.is_zero());
        assert!(ideal.contains(&Polynomial::zero(2, MonomialOrder::Lex)));
    }

    #[test]
    fn test_ideal_membership_simple() {
        let ring = ring2();
        // I = <x>
        let f = ring.from_raw(vec![(1.0, vec![1, 0])]);
        let ideal = Ideal::new(vec![f], ring.clone());
        // x*y should be in <x>
        let g = ring.from_raw(vec![(1.0, vec![1, 1])]);
        assert!(ideal.contains(&g));
        // y should NOT be in <x>
        let h = ring.from_raw(vec![(1.0, vec![0, 1])]);
        assert!(!ideal.contains(&h));
    }

    #[test]
    fn test_ideal_membership_x_squared_minus_1() {
        let ring = ring2();
        // I = <x^2 - 1>
        let f = ring.from_raw(vec![(1.0, vec![2, 0]), (-1.0, vec![0, 0])]);
        let ideal = Ideal::new(vec![f], ring.clone());
        // (x^2 - 1) * x = x^3 - x should be in the ideal
        let g = ring.from_raw(vec![(1.0, vec![3, 0]), (-1.0, vec![1, 0])]);
        assert!(ideal.contains(&g));
    }

    #[test]
    fn test_ideal_sum() {
        let ring = ring2();
        let f = ring.from_raw(vec![(1.0, vec![1, 0])]); // x
        let g = ring.from_raw(vec![(1.0, vec![0, 1])]); // y
        let i1 = Ideal::new(vec![f], ring.clone());
        let i2 = Ideal::new(vec![g], ring);
        let sum = i1.sum(&i2);
        assert_eq!(sum.num_generators(), 2);
    }

    #[test]
    fn test_ideal_product() {
        let ring = ring2();
        let f = ring.from_raw(vec![(1.0, vec![1, 0])]); // x
        let g = ring.from_raw(vec![(1.0, vec![0, 1])]); // y
        let i1 = Ideal::new(vec![f], ring.clone());
        let i2 = Ideal::new(vec![g], ring.clone());
        let prod = i1.product(&i2);
        assert_eq!(prod.num_generators(), 1);
        // product generator should be x*y
        assert!(prod.contains(&ring.from_raw(vec![(1.0, vec![1, 1])])));
    }
}
