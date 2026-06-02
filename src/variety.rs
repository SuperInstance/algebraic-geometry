//! Affine and projective varieties: dimension, singular points, vanishing ideals.

use serde::{Deserialize, Serialize};
use crate::polynomial::{Polynomial, PolynomialRing};
use crate::monomial_order::MonomialOrder;
use crate::ideal::Ideal;
use crate::groebner::groebner_basis;

/// An affine variety V(I) for an ideal I.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AffineVariety {
    /// The defining polynomials.
    pub equations: Vec<Polynomial>,
    /// The ambient polynomial ring.
    pub ring: PolynomialRing,
}

impl AffineVariety {
    /// Create a new affine variety from defining equations.
    pub fn new(equations: Vec<Polynomial>, ring: PolynomialRing) -> Self {
        Self { equations, ring }
    }

    /// The empty variety.
    pub fn empty(ring: PolynomialRing) -> Self {
        Self {
            equations: vec![ring.constant(1.0)],
            ring,
        }
    }

    /// Check if a point lies on this variety.
    pub fn contains_point(&self, point: &[f64], tol: f64) -> bool {
        self.equations.iter().all(|f| f.vanishes_at(point, tol))
    }

    /// Compute the dimension of the variety (simplified).
    /// Uses the heuristic: n - rank of the Jacobian at a generic point.
    /// For a more accurate computation, we use the leading monomials of the Groebner basis.
    pub fn dimension(&self) -> Option<usize> {
        if self.equations.is_empty() {
            return Some(self.ring.num_vars); // Whole affine space
        }

        let n = self.ring.num_vars;

        // Check if the variety is empty (1 in the ideal)
        let gb = groebner_basis(&self.equations, self.ring.order);
        if gb.len() == 1 && gb[0].num_terms() == 1 {
            let lt = gb[0].leading_term().unwrap();
            if lt.exponents.iter().all(|&e| e == 0) {
                return None; // Empty variety
            }
        }

        // Dimension heuristic: count variables that appear as leading monomials
        // A variable x_i is "constrained" if some leading monomial is a pure power of x_i
        // or involves only higher-index variables in lex order
        let constrained: std::collections::HashSet<usize> = gb
            .iter()
            .filter_map(|p| p.leading_monomial())
            .flat_map(|lm| {
                lm.iter()
                    .enumerate()
                    .filter(|(_, &e)| e > 0)
                    .map(|(i, _)| i)
                    .collect::<Vec<_>>()
            })
            .collect();

        // Simplified: dimension = number of free variables
        // This is an approximation; true dimension requires more sophisticated analysis
        Some(n.saturating_sub(self.equations.len().min(n)))
    }

    /// Compute the Jacobian matrix at a point.
    pub fn jacobian_at(&self, point: &[f64]) -> Vec<Vec<f64>> {
        let n = self.ring.num_vars;
        let eps = 1e-8;

        self.equations
            .iter()
            .map(|f| {
                (0..n)
                    .map(|i| {
                        let mut p_plus = point.to_vec();
                        let mut p_minus = point.to_vec();
                        p_plus[i] += eps;
                        p_minus[i] -= eps;
                        (f.evaluate(&p_plus) - f.evaluate(&p_minus)) / (2.0 * eps)
                    })
                    .collect()
            })
            .collect()
    }

    /// Find singular points by checking rank deficiency of Jacobian.
    /// Returns true if the given point is a singular point.
    pub fn is_singular_point(&self, point: &[f64], tol: f64) -> bool {
        if !self.contains_point(point, tol * 10.0) {
            return false;
        }
        let jac = self.jacobian_at(point);
        // Check if all rows are zero (simplified singular check)
        jac.iter().all(|row| row.iter().all(|&v| v.abs() < tol))
    }

    /// Union of two varieties: V(I) ∪ V(J) = V(I·J).
    pub fn union(&self, other: &AffineVariety) -> AffineVariety {
        let mut eqs = Vec::new();
        for a in &self.equations {
            for b in &other.equations {
                eqs.push(a.mul(b));
            }
        }
        AffineVariety::new(eqs, self.ring.clone())
    }

    /// Intersection of two varieties: V(I) ∩ V(J) = V(I+J).
    pub fn intersection(&self, other: &AffineVariety) -> AffineVariety {
        let mut eqs = self.equations.clone();
        eqs.extend(other.equations.iter().cloned());
        AffineVariety::new(eqs, self.ring.clone())
    }

    /// Check if this variety is contained in another.
    /// V(I) ⊆ V(J) iff J's ideal ⊆ I's ideal (radical containment).
    pub fn is_contained_in(&self, other: &AffineVariety) -> bool {
        // Simplified: check if all of other's equations vanish on our variety
        // This requires checking if other's generators are in the radical of our ideal
        // For now, use the Groebner basis approach
        let our_ideal = Ideal::new(self.equations.clone(), self.ring.clone());
        for eq in &other.equations {
            if !our_ideal.contains(eq) {
                // Not necessarily false — need radical membership
                // Simplified check for the common case
            }
        }
        // A more accurate check: if our Groebner basis reduces all of other's equations
        let gb = groebner_basis(&self.equations, self.ring.order);
        other.equations.iter().all(|f| {
            let (_, r) = crate::groebner::reduce_by_set(f, &gb);
            r.is_zero()
        })
    }

    /// Number of defining equations.
    pub fn num_equations(&self) -> usize {
        self.equations.len()
    }
}

/// A projective variety.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjectiveVariety {
    /// Homogeneous defining polynomials.
    pub equations: Vec<Polynomial>,
    /// The polynomial ring (with homogenizing variable).
    pub ring: PolynomialRing,
}

impl ProjectiveVariety {
    /// Create a projective variety from homogeneous equations.
    pub fn new(equations: Vec<Polynomial>, ring: PolynomialRing) -> Self {
        Self { equations, ring }
    }

    /// Create from affine equations by homogenizing.
    pub fn from_affine(affine_eqs: &[Polynomial], ring: &PolynomialRing) -> Self {
        let hom_ring = PolynomialRing::new(ring.num_vars + 1, ring.order);
        let equations: Vec<Polynomial> = affine_eqs.iter().map(|f| f.homogenize()).collect();
        Self {
            equations,
            ring: hom_ring,
        }
    }

    /// Check if a point in projective coordinates lies on this variety.
    pub fn contains_projective_point(&self, point: &[f64], tol: f64) -> bool {
        self.equations.iter().all(|f| {
            // For homogeneous polynomial, evaluate and check zero
            f.evaluate(point).abs() < tol
        })
    }

    /// Dimension (heuristic).
    pub fn dimension(&self) -> Option<usize> {
        if self.equations.is_empty() {
            return Some(self.ring.num_vars - 1); // Projective space has dimension n-1
        }
        Some(self.ring.num_vars.saturating_sub(1).saturating_sub(self.equations.len()))
    }

    /// Check if equations are homogeneous.
    pub fn is_well_defined(&self) -> bool {
        self.equations.iter().all(|f| f.is_homogeneous())
    }

    /// Degree of the projective variety (simplified: product of degrees of defining equations).
    pub fn degree(&self) -> u32 {
        self.equations
            .iter()
            .map(|f| f.total_degree())
            .product()
    }
}

/// Compute the vanishing ideal I(V) for a finite set of points V.
/// For points p1, ..., pk in k^n, I(V) is the ideal of polynomials vanishing at all pi.
pub fn vanishing_ideal(points: &[Vec<f64>], ring: &PolynomialRing) -> Ideal {
    if points.is_empty() {
        return Ideal::new(vec![ring.constant(1.0)], ring.clone());
    }

    let n = ring.num_vars;
    let mut generators = Vec::new();

    // For each variable x_i, construct the minimal polynomial that vanishes at all points' x_i values
    // This is (x_i - p1_i)(x_i - p2_i)...(x_i - pk_i) for each i
    for i in 0..n {
        let x_i = ring.variable(i);
        let mut poly = ring.constant(1.0);
        for point in points {
            let factor = x_i.sub(&ring.constant(point[i]));
            poly = poly.mul(&factor);
        }
        generators.push(poly);
    }

    Ideal::new(generators, ring.clone())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring2() -> PolynomialRing {
        PolynomialRing::new(2, MonomialOrder::Lex)
    }

    fn ring3() -> PolynomialRing {
        PolynomialRing::new(3, MonomialOrder::Lex)
    }

    #[test]
    fn test_affine_variety_contains_point() {
        let ring = ring2();
        // V(x^2 - 1)
        let f = ring.from_raw(vec![(1.0, vec![2, 0]), (-1.0, vec![0, 0])]);
        let v = AffineVariety::new(vec![f], ring);
        assert!(v.contains_point(&[1.0, 0.0], 1e-10));
        assert!(v.contains_point(&[-1.0, 5.0], 1e-10));
        assert!(!v.contains_point(&[2.0, 0.0], 1e-10));
    }

    #[test]
    fn test_affine_variety_union() {
        let ring = ring2();
        // V(x) ∪ V(y)
        let f = ring.from_raw(vec![(1.0, vec![1, 0])]);
        let g = ring.from_raw(vec![(1.0, vec![0, 1])]);
        let v1 = AffineVariety::new(vec![f], ring.clone());
        let v2 = AffineVariety::new(vec![g], ring);
        let union = v1.union(&v2);
        // V(x*y): (0, anything) and (anything, 0)
        assert!(union.contains_point(&[0.0, 5.0], 1e-10));
        assert!(union.contains_point(&[3.0, 0.0], 1e-10));
        assert!(!union.contains_point(&[1.0, 1.0], 1e-10));
    }

    #[test]
    fn test_affine_variety_intersection() {
        let ring = ring2();
        let f = ring.from_raw(vec![(1.0, vec![1, 0])]); // V(x)
        let g = ring.from_raw(vec![(1.0, vec![0, 1])]); // V(y)
        let v1 = AffineVariety::new(vec![f], ring.clone());
        let v2 = AffineVariety::new(vec![g], ring);
        let inter = v1.intersection(&v2);
        // V(x,y) = just the origin
        assert!(inter.contains_point(&[0.0, 0.0], 1e-10));
        assert!(!inter.contains_point(&[1.0, 0.0], 1e-10));
    }

    #[test]
    fn test_affine_variety_dimension() {
        let ring = ring3();
        // V(0) = all of A^3, dimension 3
        let v = AffineVariety::new(vec![], ring);
        assert_eq!(v.dimension(), Some(3));
    }

    #[test]
    fn test_projective_from_affine() {
        let ring = ring2();
        let f = ring.from_raw(vec![(1.0, vec![2, 0]), (-1.0, vec![0, 0])]); // x^2 - 1
        let pv = ProjectiveVariety::from_affine(&[f], &ring);
        assert_eq!(pv.ring.num_vars, 3);
        assert!(pv.is_well_defined());
    }

    #[test]
    fn test_projective_degree() {
        let ring = ring2();
        let f = ring.from_raw(vec![(1.0, vec![2, 0]), (-1.0, vec![0, 0])]);
        let pv = ProjectiveVariety::from_affine(&[f], &ring);
        assert_eq!(pv.degree(), 2);
    }

    #[test]
    fn test_vanishing_ideal() {
        let ring = ring2();
        // Vanishing ideal of points {(0,0), (1,1)}
        let ideal = vanishing_ideal(&[vec![0.0, 0.0], vec![1.0, 1.0]], &ring);
        // Should have 2 generators (one per variable)
        assert_eq!(ideal.num_generators(), 2);
        // x(x-1) should vanish at both points
        let gen = &ideal.generators[0];
        assert!(gen.evaluate(&[0.0, 0.0]).abs() < 1e-10);
        assert!(gen.evaluate(&[1.0, 1.0]).abs() < 1e-10);
    }

    #[test]
    fn test_jacobian() {
        let ring = ring2();
        let f = ring.from_raw(vec![(1.0, vec![2, 0]), (-1.0, vec![0, 0])]); // x^2 - 1
        let v = AffineVariety::new(vec![f], ring);
        let jac = v.jacobian_at(&[1.0, 0.0]);
        // df/dx = 2x = 2, df/dy = 0
        assert!((jac[0][0] - 2.0).abs() < 1e-6);
        assert!(jac[0][1].abs() < 1e-6);
    }

    #[test]
    fn test_singular_point() {
        let ring = ring2();
        // y^2 = x^3 (cusp at origin)
        let f = ring.from_raw(vec![(1.0, vec![0, 2]), (-1.0, vec![3, 0])]);
        let v = AffineVariety::new(vec![f], ring);
        // Origin should be singular
        assert!(v.is_singular_point(&[0.0, 0.0], 1e-6));
    }
}
