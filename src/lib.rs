//! # algebraic-geometry
//!
//! Computational algebraic geometry: polynomial rings, ideals, varieties,
//! Groebner bases (Buchberger's algorithm), elimination theory, Bezout's theorem,
//! and polynomial constraint solving.

pub mod polynomial;
pub mod monomial_order;
pub mod ideal;
pub mod groebner;
pub mod variety;
pub mod elimination;
pub mod bezout;
pub mod constraint_solving;

pub use polynomial::{Polynomial, PolynomialRing};
pub use monomial_order::MonomialOrder;
pub use ideal::Ideal;
pub use groebner::groebner_basis;
pub use variety::{AffineVariety, ProjectiveVariety};
pub use elimination::elimination_ideal;
pub use bezout::bezout_count;
pub use constraint_solving::solve_constraints;
