//! Multivariate polynomials over a field (rationals via f64 approximations for simplicity).

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;

use crate::monomial_order::MonomialOrder;

/// A monomial represented by its exponent vector and a coefficient.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Term {
    pub coeff: f64,
    pub exponents: Vec<u32>,
}

impl Term {
    pub fn new(coeff: f64, exponents: Vec<u32>) -> Self {
        Self { coeff, exponents }
    }

    /// Total degree of this term.
    pub fn total_degree(&self) -> u32 {
        self.exponents.iter().sum()
    }

    /// Multiply two terms (add exponents, multiply coefficients).
    pub fn multiply(&self, other: &Term) -> Term {
        let exponents: Vec<u32> = self
            .exponents
            .iter()
            .zip(other.exponents.iter())
            .map(|(a, b)| a + b)
            .collect();
        Term::new(self.coeff * other.coeff, exponents)
    }

    /// Check if this term divides other (exponent-wise).
    pub fn divides(&self, other: &Term) -> bool {
        self.exponents
            .iter()
            .zip(other.exponents.iter())
            .all(|(a, b)| a <= b)
    }

    /// Compute other / self (assuming self divides other).
    pub fn quotient(&self, other: &Term) -> Term {
        let exponents: Vec<u32> = other
            .exponents
            .iter()
            .zip(self.exponents.iter())
            .map(|(a, b)| a - b)
            .collect();
        Term::new(other.coeff / self.coeff, exponents)
    }
}

/// A multivariate polynomial.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Polynomial {
    /// Terms stored as map from exponent vector (as tuple-keyed) to coefficient.
    /// We store as Vec<Term> for simplicity.
    pub terms: Vec<Term>,
    /// Number of variables.
    pub num_vars: usize,
    /// The monomial ordering used for this polynomial.
    pub order: MonomialOrder,
}

impl Polynomial {
    /// Create a new zero polynomial with given number of variables and order.
    pub fn zero(num_vars: usize, order: MonomialOrder) -> Self {
        Self {
            terms: Vec::new(),
            num_vars,
            order,
        }
    }

    /// Create a constant polynomial.
    pub fn constant(c: f64, num_vars: usize, order: MonomialOrder) -> Self {
        if c.abs() < 1e-12 {
            return Self::zero(num_vars, order);
        }
        Self {
            terms: vec![Term::new(c, vec![0; num_vars])],
            num_vars,
            order,
        }
    }

    /// Create a polynomial from a single term.
    pub fn from_term(term: Term, num_vars: usize, order: MonomialOrder) -> Self {
        let mut p = Self::zero(num_vars, order);
        if term.coeff.abs() > 1e-12 {
            p.terms.push(term);
        }
        p
    }

    /// Create a polynomial representing a single variable x_i.
    pub fn variable(var_index: usize, num_vars: usize, order: MonomialOrder) -> Self {
        let mut exp = vec![0u32; num_vars];
        exp[var_index] = 1;
        Self {
            terms: vec![Term::new(1.0, exp)],
            num_vars,
            order,
        }
    }

    /// Is this the zero polynomial?
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Leading term (under the polynomial's monomial order).
    pub fn leading_term(&self) -> Option<&Term> {
        if self.terms.is_empty() {
            return None;
        }
        let mut best = &self.terms[0];
        for term in &self.terms[1..] {
            if self.order.compare(&term.exponents, &best.exponents) == std::cmp::Ordering::Greater {
                best = term;
            }
        }
        Some(best)
    }

    /// Leading monomial (exponent vector of leading term).
    pub fn leading_monomial(&self) -> Option<Vec<u32>> {
        self.leading_term().map(|t| t.exponents.clone())
    }

    /// Leading coefficient.
    pub fn leading_coefficient(&self) -> Option<f64> {
        self.leading_term().map(|t| t.coeff)
    }

    /// Total degree of the polynomial.
    pub fn total_degree(&self) -> u32 {
        self.terms.iter().map(|t| t.total_degree()).max().unwrap_or(0)
    }

    /// Collect like terms and remove near-zero coefficients.
    pub fn simplify(&mut self) {
        let mut map: HashMap<Vec<u32>, f64> = HashMap::new();
        for term in self.terms.drain(..) {
            *map.entry(term.exponents).or_insert(0.0) += term.coeff;
        }
        self.terms = map
            .into_iter()
            .filter(|(_, c)| c.abs() > 1e-12)
            .map(|(exponents, coeff)| Term::new(coeff, exponents))
            .collect();
    }

    /// Add two polynomials.
    pub fn add(&self, other: &Polynomial) -> Polynomial {
        let mut result = self.clone();
        result.terms.extend(other.terms.iter().cloned());
        result.simplify();
        result
    }

    /// Subtract polynomials.
    pub fn sub(&self, other: &Polynomial) -> Polynomial {
        let neg: Polynomial = Polynomial {
            terms: other
                .terms
                .iter()
                .map(|t| Term::new(-t.coeff, t.exponents.clone()))
                .collect(),
            num_vars: other.num_vars,
            order: other.order,
        };
        self.add(&neg)
    }

    /// Multiply two polynomials.
    pub fn mul(&self, other: &Polynomial) -> Polynomial {
        let mut result = Polynomial::zero(self.num_vars, self.order);
        for a in &self.terms {
            for b in &other.terms {
                result.terms.push(a.multiply(b));
            }
        }
        result.simplify();
        result
    }

    /// Scalar multiplication.
    pub fn scalar_mul(&self, c: f64) -> Polynomial {
        Polynomial {
            terms: self
                .terms
                .iter()
                .map(|t| Term::new(t.coeff * c, t.exponents.clone()))
                .filter(|t| t.coeff.abs() > 1e-12)
                .collect(),
            num_vars: self.num_vars,
            order: self.order,
        }
    }

    /// Evaluate the polynomial at a given point.
    pub fn evaluate(&self, point: &[f64]) -> f64 {
        self.terms
            .iter()
            .map(|term| {
                let mut val = term.coeff;
                for (i, &e) in term.exponents.iter().enumerate() {
                    val *= point[i].powi(e as i32);
                }
                val
            })
            .sum()
    }

    /// Check if a point is approximately a root.
    pub fn vanishes_at(&self, point: &[f64], tol: f64) -> bool {
        self.evaluate(point).abs() < tol
    }

    /// Divide polynomial by a single term, returning quotient or None.
    pub fn divide_by_term(&self, divisor: &Term) -> Option<Polynomial> {
        let mut result = Polynomial::zero(self.num_vars, self.order);
        for term in &self.terms {
            if divisor.divides(term) {
                result.terms.push(divisor.quotient(term));
            } else {
                return None;
            }
        }
        result.simplify();
        Some(result)
    }

    /// Check if the polynomial is homogeneous.
    pub fn is_homogeneous(&self) -> bool {
        if self.terms.is_empty() {
            return true;
        }
        let deg = self.terms[0].total_degree();
        self.terms.iter().all(|t| t.total_degree() == deg)
    }

    /// Homogenize a polynomial by adding a new variable.
    pub fn homogenize(&self) -> Polynomial {
        let deg = self.total_degree();
        let mut result = Polynomial::zero(self.num_vars + 1, self.order);
        for term in &self.terms {
            let mut exp = term.exponents.clone();
            exp.push(deg - term.total_degree());
            result.terms.push(Term::new(term.coeff, exp));
        }
        result
    }

    /// Number of terms.
    pub fn num_terms(&self) -> usize {
        self.terms.len()
    }
}

impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.terms.is_empty() {
            return write!(f, "0");
        }
        let var_names: Vec<char> = "xyzwuvstrabcdefgh".chars().collect();
        let mut parts = Vec::new();
        for term in &self.terms {
            let mut s = String::new();
            if term.coeff.abs() > 1e-12 {
                let has_vars = term.exponents.iter().any(|&e| e > 0);
                if !has_vars {
                    s.push_str(&format!("{}", term.coeff));
                } else if (term.coeff - 1.0).abs() < 1e-12 {
                    // coeff is 1, skip
                } else if (term.coeff + 1.0).abs() < 1e-12 {
                    s.push('-');
                } else {
                    s.push_str(&format!("{}*", term.coeff));
                }
                for (i, &e) in term.exponents.iter().enumerate() {
                    if e > 0 {
                        if i < var_names.len() {
                            s.push(var_names[i]);
                        } else {
                            s.push_str(&format!("x{}", i));
                        }
                        if e > 1 {
                            s.push_str(&format!("^{}", e));
                        }
                    }
                }
                if !s.is_empty() {
                    parts.push(s);
                }
            }
        }
        if parts.is_empty() {
            write!(f, "0")
        } else {
            write!(f, "{}", parts.join(" + ").replace("+ -", "- "))
        }
    }
}

/// A polynomial ring with specified number of variables and monomial order.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolynomialRing {
    pub num_vars: usize,
    pub var_names: Vec<String>,
    pub order: MonomialOrder,
}

impl PolynomialRing {
    pub fn new(num_vars: usize, order: MonomialOrder) -> Self {
        let default_names: Vec<String> = "xyzwuvstrabcdefgh"
            .chars()
            .take(num_vars)
            .map(|c| c.to_string())
            .collect();
        let var_names = if num_vars <= default_names.len() {
            default_names
        } else {
            (0..num_vars).map(|i| format!("x{}", i)).collect()
        };
        Self {
            num_vars,
            var_names,
            order,
        }
    }

    pub fn zero(&self) -> Polynomial {
        Polynomial::zero(self.num_vars, self.order)
    }

    pub fn constant(&self, c: f64) -> Polynomial {
        Polynomial::constant(c, self.num_vars, self.order)
    }

    pub fn variable(&self, idx: usize) -> Polynomial {
        Polynomial::variable(idx, self.num_vars, self.order)
    }

    /// Create polynomial from (coeff, exponents) pairs.
    pub fn from_raw(&self, terms: Vec<(f64, Vec<u32>)>) -> Polynomial {
        let mut p = Polynomial::zero(self.num_vars, self.order);
        for (coeff, exp) in terms {
            if coeff.abs() > 1e-12 {
                p.terms.push(Term::new(coeff, exp));
            }
        }
        p.simplify();
        p
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring2() -> PolynomialRing {
        PolynomialRing::new(2, MonomialOrder::Lex)
    }

    fn ring3() -> PolynomialRing {
        PolynomialRing::new(3, MonomialOrder::GrLex)
    }

    #[test]
    fn test_zero_polynomial() {
        let p = Polynomial::zero(2, MonomialOrder::Lex);
        assert!(p.is_zero());
        assert_eq!(p.total_degree(), 0);
        assert!(p.leading_term().is_none());
    }

    #[test]
    fn test_constant_polynomial() {
        let p = Polynomial::constant(5.0, 2, MonomialOrder::Lex);
        assert!(!p.is_zero());
        assert_eq!(p.evaluate(&[0.0, 0.0]), 5.0);
    }

    #[test]
    fn test_variable_polynomial() {
        let p = Polynomial::variable(0, 2, MonomialOrder::Lex); // x
        assert_eq!(p.evaluate(&[3.0, 0.0]), 3.0);
        assert_eq!(p.evaluate(&[0.0, 5.0]), 0.0);
    }

    #[test]
    fn test_polynomial_add() {
        let ring = ring2();
        let p1 = ring.from_raw(vec![(1.0, vec![2, 0]), (1.0, vec![0, 1])]); // x^2 + y
        let p2 = ring.from_raw(vec![(1.0, vec![2, 0]), (2.0, vec![1, 0])]); // x^2 + 2x
        let sum = p1.add(&p2);
        assert_eq!(sum.terms.len(), 3);
        // 2x^2 + y + 2x
        assert_eq!(sum.evaluate(&[1.0, 1.0]), 5.0);
    }

    #[test]
    fn test_polynomial_sub() {
        let ring = ring2();
        let p1 = ring.from_raw(vec![(3.0, vec![1, 0])]); // 3x
        let p2 = ring.from_raw(vec![(1.0, vec![1, 0])]); // x
        let diff = p1.sub(&p2);
        assert_eq!(diff.evaluate(&[5.0, 0.0]), 10.0);
    }

    #[test]
    fn test_polynomial_mul() {
        let ring = ring2();
        // (x + 1)(x - 1) = x^2 - 1
        let p1 = ring.from_raw(vec![(1.0, vec![1, 0]), (1.0, vec![0, 0])]);
        let p2 = ring.from_raw(vec![(1.0, vec![1, 0]), (-1.0, vec![0, 0])]);
        let prod = p1.mul(&p2);
        assert!((prod.evaluate(&[2.0, 0.0]) - 3.0).abs() < 1e-10);
        assert!((prod.evaluate(&[1.0, 0.0])).abs() < 1e-10);
    }

    #[test]
    fn test_scalar_mul() {
        let ring = ring2();
        let p = ring.from_raw(vec![(1.0, vec![1, 0]), (2.0, vec![0, 1])]);
        let scaled = p.scalar_mul(3.0);
        assert_eq!(scaled.evaluate(&[1.0, 1.0]), 9.0);
    }

    #[test]
    fn test_leading_term_lex() {
        let ring = ring2();
        let p = ring.from_raw(vec![
            (1.0, vec![0, 2]),  // y^2
            (1.0, vec![1, 0]),  // x
            (1.0, vec![2, 0]),  // x^2
        ]);
        let lt = p.leading_term().unwrap();
        assert_eq!(lt.exponents, vec![2, 0]); // x^2 is leading under lex
    }

    #[test]
    fn test_leading_term_grlex() {
        let ring = PolynomialRing::new(2, MonomialOrder::GrLex);
        let p = ring.from_raw(vec![
            (1.0, vec![2, 0]),  // x^2, deg 2
            (1.0, vec![1, 1]),  // x*y, deg 2
            (1.0, vec![0, 3]),  // y^3, deg 3
        ]);
        let lt = p.leading_term().unwrap();
        assert_eq!(lt.exponents, vec![0, 3]); // y^3 has highest degree
    }

    #[test]
    fn test_vanishes_at() {
        let ring = ring2();
        let p = ring.from_raw(vec![(1.0, vec![2, 0]), (-1.0, vec![0, 0])]); // x^2 - 1
        assert!(p.vanishes_at(&[1.0, 0.0], 1e-10));
        assert!(p.vanishes_at(&[-1.0, 0.0], 1e-10));
        assert!(!p.vanishes_at(&[2.0, 0.0], 1e-10));
    }

    #[test]
    fn test_term_divides() {
        let t1 = Term::new(1.0, vec![1, 0]); // x
        let t2 = Term::new(1.0, vec![2, 1]); // x^2*y
        assert!(t1.divides(&t2));
        assert!(!t2.divides(&t1));
    }

    #[test]
    fn test_ring_from_raw() {
        let ring = ring3();
        // x + 2y + z^3
        let p = ring.from_raw(vec![
            (1.0, vec![1, 0, 0]),
            (2.0, vec![0, 1, 0]),
            (1.0, vec![0, 0, 3]),
        ]);
        assert_eq!(p.num_terms(), 3);
        assert_eq!(p.evaluate(&[1.0, 1.0, 1.0]), 4.0);
    }

    #[test]
    fn test_homogeneous() {
        let ring = ring2();
        let p = ring.from_raw(vec![
            (1.0, vec![2, 0]),
            (2.0, vec![1, 1]),
            (1.0, vec![0, 2]),
        ]);
        assert!(p.is_homogeneous());
    }

    #[test]
    fn test_not_homogeneous() {
        let ring = ring2();
        let p = ring.from_raw(vec![(1.0, vec![2, 0]), (1.0, vec![0, 0])]);
        assert!(!p.is_homogeneous());
    }

    #[test]
    fn test_homogenize() {
        let ring = ring2();
        // x^2 + y -> x^2 + y*z (homogenized with z)
        let p = ring.from_raw(vec![(1.0, vec![2, 0]), (1.0, vec![0, 1])]);
        let h = p.homogenize();
        assert!(h.is_homogeneous());
        assert_eq!(h.num_vars, 3);
    }

    #[test]
    fn test_display() {
        let ring = ring2();
        let p = ring.from_raw(vec![(1.0, vec![2, 0]), (1.0, vec![0, 1])]);
        let s = format!("{}", p);
        assert!(!s.is_empty());
    }
}
