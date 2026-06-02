//! Application: polynomial constraint solving — polynomial systems from polynomial constraints.
//!
//! Models polynomial systems as polynomial constraint systems and solves them using
//! Groebner bases and elimination theory.

use crate::monomial_order::MonomialOrder;
use crate::polynomial::{Polynomial, PolynomialRing};
use crate::groebner::groebner_basis;
use crate::elimination::solve_elimination;
use crate::ideal::Ideal;
use serde::{Deserialize, Serialize};

/// An polynomial constraint expressed as a polynomial equation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolynomialConstraint {
    /// Human-readable description of the constraint.
    pub description: String,
    /// The polynomial representing the constraint.
    pub polynomial: Polynomial,
    /// Constraint type.
    pub constraint_type: ConstraintKind,
}

/// Types of constraints.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ConstraintKind {
    /// Equality constraint: f(x) = 0.
    Equality,
    /// Inequality constraint: f(x) > 0 (approximated).
    Inequality,
    /// Resource constraint: sum of resources <= limit.
    ResourceBound,
    /// Custom constraint: polynomial systems must satisfy some property.
    Custom,
}

/// A system of polynomial constraints.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolynomialConstraintSystem {
    pub constraints: Vec<PolynomialConstraint>,
    pub ring: PolynomialRing,
    /// Variable names (semantic labels).
    pub variable_labels: Vec<String>,
}

impl PolynomialConstraintSystem {
    /// Create a new constraint system.
    pub fn new(num_vars: usize, variable_labels: Vec<String>) -> Self {
        Self {
            constraints: Vec::new(),
            ring: PolynomialRing::new(num_vars, MonomialOrder::Lex),
            variable_labels,
        }
    }

    /// Add an equality constraint.
    pub fn add_equality(&mut self, description: &str, polynomial: Polynomial) {
        self.constraints.push(PolynomialConstraint {
            description: description.to_string(),
            polynomial,
            constraint_type: ConstraintKind::Equality,
        });
    }

    /// Add a resource constraint (sum_i coeff_i * x_i <= limit).
    pub fn add_resource_constraint(&mut self, coeffs: &[f64], limit: f64, description: &str) {
        let terms: Vec<(f64, Vec<u32>)> = coeffs
            .iter()
            .enumerate()
            .filter(|(_, &c)| c.abs() > 1e-12)
            .map(|(i, &c)| {
                let mut exp = vec![0u32; self.ring.num_vars];
                exp[i] = 1;
                (c, exp)
            })
            .collect();
        let poly = self.ring.from_raw(terms);
        let constraint_poly = poly.sub(&self.ring.constant(limit));
        self.constraints.push(PolynomialConstraint {
            description: description.to_string(),
            polynomial: constraint_poly,
            constraint_type: ConstraintKind::ResourceBound,
        });
    }

    /// Get all equality constraints as polynomials.
    pub fn equality_polynomials(&self) -> Vec<Polynomial> {
        self.constraints
            .iter()
            .filter(|c| c.constraint_type == ConstraintKind::Equality)
            .map(|c| c.polynomial.clone())
            .collect()
    }

    /// Check if a point satisfies all constraints.
    pub fn satisfies(&self, point: &[f64], tol: f64) -> bool {
        self.constraints.iter().all(|c| {
            let val = c.polynomial.evaluate(point);
            match c.constraint_type {
                ConstraintKind::Equality | ConstraintKind::ResourceBound => val.abs() < tol,
                ConstraintKind::Inequality => val > -tol,
                ConstraintKind::Custom => val.abs() < tol,
            }
        })
    }
}

/// Result of solving a constraint system.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConstraintSolution {
    /// Solutions found (variable values).
    pub solutions: Vec<Vec<f64>>,
    /// Whether the system is consistent.
    pub consistent: bool,
    /// Number of solutions found.
    pub num_solutions: usize,
    /// The Groebner basis used.
    pub groebner_basis: Vec<Polynomial>,
}

/// Solve a constraint system using Groebner bases.
pub fn solve_constraints(system: &PolynomialConstraintSystem) -> ConstraintSolution {
    let polys = system.equality_polynomials();

    if polys.is_empty() {
        return ConstraintSolution {
            solutions: vec![],
            consistent: true,
            num_solutions: 0,
            groebner_basis: vec![],
        };
    }

    let gb = groebner_basis(&polys, MonomialOrder::Lex);

    // Check for inconsistency: if 1 is in the Groebner basis
    let inconsistent = gb.iter().any(|p| {
        p.num_terms() == 1 && p.terms[0].exponents.iter().all(|&e| e == 0) && p.terms[0].coeff.abs() > 0.5
    });

    if inconsistent {
        return ConstraintSolution {
            solutions: vec![],
            consistent: false,
            num_solutions: 0,
            groebner_basis: gb,
        };
    }

    let solutions = solve_elimination(&polys, system.ring.num_vars);

    ConstraintSolution {
        num_solutions: solutions.len(),
        consistent: true,
        solutions,
        groebner_basis: gb,
    }
}

/// Model a simple resource allocation as polynomial constraints.
/// Each variable x_i represents an allocation.
/// Constraints: sum = total, each x_i >= 0 (approximated).
pub fn resource_allocation(num_vars: usize, total: f64, min_allocations: &[f64]) -> PolynomialConstraintSystem {
    let mut system = PolynomialConstraintSystem::new(
        num_vars,
        (0..num_vars).map(|i| format!("var_{}", i)).collect(),
    );

    // Constraint: sum of allocations = total
    let coeffs: Vec<f64> = vec![1.0; num_vars];
    system.add_resource_constraint(&coeffs, total, "total resource allocation");

    // Add minimum value constraints: x_i >= min_i => x_i - min_i = 0 (simplified)
    for (i, &min_alloc) in min_allocations.iter().enumerate() {
        if min_alloc.abs() > 1e-12 {
            let x_i = system.ring.variable(i);
            let constraint = x_i.sub(&system.ring.constant(min_alloc));
            system.add_equality(
                &format!("var_{} minimum value", i),
                constraint,
            );
        }
    }

    system
}

/// Model behavioral constraints: var_i produces output proportional to input^power.
pub fn production_constraints(
    num_vars: usize,
    powers: &[f64],
    total_input: f64,
    total_output: f64,
) -> PolynomialConstraintSystem {
    let mut system = PolynomialConstraintSystem::new(
        num_vars,
        (0..num_vars).map(|i| format!("input_{}", i)).collect(),
    );

    // Total input constraint
    let coeffs: Vec<f64> = vec![1.0; num_vars];
    system.add_resource_constraint(&coeffs, total_input, "total input");

    // Production: sum of x_i^power_i = total_output (approximated as linear for simplicity)
    // For integer powers, we can represent exactly
    let terms: Vec<(f64, Vec<u32>)> = powers
        .iter()
        .enumerate()
        .map(|(i, &p)| {
            let mut exp = vec![0u32; num_vars];
            exp[i] = p.round() as u32;
            (1.0, exp)
        })
        .collect();
    let prod_poly = system.ring.from_raw(terms);
    let constraint = prod_poly.sub(&system.ring.constant(total_output));
    system.add_equality("total output constraint", constraint);

    system
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_var_system_creation() {
        let system = PolynomialConstraintSystem::new(3, vec!["x".into(), "y".into(), "z".into()]);
        assert_eq!(system.constraints.len(), 0);
        assert_eq!(system.ring.num_vars, 3);
    }

    #[test]
    fn test_var_equality_constraint() {
        let mut system = PolynomialConstraintSystem::new(2, vec!["x".into(), "y".into()]);
        // x + y = 10
        let poly = system.ring.from_raw(vec![
            (1.0, vec![1, 0]),
            (1.0, vec![0, 1]),
        ]);
        system.add_equality("x + y = 10", poly.sub(&system.ring.constant(10.0)));
        assert_eq!(system.constraints.len(), 1);
    }

    #[test]
    fn test_solve_simple_system() {
        let mut system = PolynomialConstraintSystem::new(2, vec!["x".into(), "y".into()]);
        // x + y = 0
        let f1 = system.ring.from_raw(vec![(1.0, vec![1, 0]), (1.0, vec![0, 1])]);
        system.add_equality("x + y = 0", f1);
        // x - y = 0
        let f2 = system.ring.from_raw(vec![(1.0, vec![1, 0]), (-1.0, vec![0, 1])]);
        system.add_equality("x - y = 0", f2);

        let result = solve_constraints(&system);
        assert!(result.consistent);
        assert!(!result.groebner_basis.is_empty());
    }

    #[test]
    fn test_inconsistent_system() {
        let mut system = PolynomialConstraintSystem::new(2, vec!["x".into(), "y".into()]);
        // x = 0
        system.add_equality("x = 0", system.ring.variable(0));
        // x = 1
        let x_minus_1 = system.ring.variable(0).sub(&system.ring.constant(1.0));
        system.add_equality("x = 1", x_minus_1);

        let result = solve_constraints(&system);
        assert!(!result.consistent);
    }

    #[test]
    fn test_resource_allocation() {
        let system = resource_allocation(3, 10.0, &[2.0, 3.0, 5.0]);
        assert!(system.constraints.len() >= 2);
    }

    #[test]
    fn test_satisfies_point() {
        let mut system = PolynomialConstraintSystem::new(2, vec!["x".into(), "y".into()]);
        // x + y = 10
        let poly = system
            .ring
            .from_raw(vec![(1.0, vec![1, 0]), (1.0, vec![0, 1])])
            .sub(&system.ring.constant(10.0));
        system.add_equality("sum = 10", poly);

        assert!(system.satisfies(&[5.0, 5.0], 1e-6));
        assert!(!system.satisfies(&[3.0, 3.0], 1e-6));
    }

    #[test]
    fn test_production_constraints() {
        let system = production_constraints(2, &[2.0, 2.0], 10.0, 50.0);
        assert!(system.constraints.len() >= 2);
    }

    #[test]
    fn test_groebner_basis_computed() {
        let mut system = PolynomialConstraintSystem::new(2, vec!["x".into(), "y".into()]);
        let f1 = system.ring.from_raw(vec![(1.0, vec![1, 0]), (1.0, vec![0, 1])]);
        system.add_equality("x + y", f1);

        let result = solve_constraints(&system);
        assert!(!result.groebner_basis.is_empty());
    }
}
