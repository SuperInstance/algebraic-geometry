//! Elimination theory: elimination ideals, extension theorem.

use crate::monomial_order::MonomialOrder;
use crate::polynomial::Polynomial;
use crate::groebner::groebner_basis;

/// Compute the l-th elimination ideal I_l = I ∩ k[x_{l+1}, ..., x_n].
/// Uses lex order Groebner basis and selects polynomials not involving x_0, ..., x_{l-1}.
pub fn elimination_ideal(generators: &[Polynomial], l: usize) -> Vec<Polynomial> {
    if generators.is_empty() {
        return vec![];
    }
    let num_vars = generators[0].num_vars;
    if l >= num_vars {
        return vec![];
    }

    // Compute Groebner basis with lex order (required for elimination)
    let gb = groebner_basis(generators, MonomialOrder::Lex);

    // Select elements that don't involve x_0, ..., x_{l-1}
    gb.into_iter()
        .filter(|p| {
            p.terms.iter().all(|term| {
                term.exponents[..l].iter().all(|&e| e == 0)
            })
        })
        .collect()
}

/// Extension theorem check: given a partial solution (x_{l+1}, ..., x_n),
/// check if it can be extended to a full solution.
/// For the simplified version, we check if the partial solution satisfies the elimination ideal.
pub fn can_extend(
    generators: &[Polynomial],
    l: usize,
    partial_point: &[f64],
) -> bool {
    if generators.is_empty() {
        return true;
    }
    let num_vars = generators[0].num_vars;
    let elim = elimination_ideal(generators, l);

    // Check that the partial point satisfies all elimination ideal polynomials
    elim.iter().all(|p| {
        let mut full_point = vec![0.0; num_vars];
        full_point[l..].copy_from_slice(partial_point);
        p.evaluate(&full_point).abs() < 1e-8
    })
}

/// Solve a polynomial system using elimination.
/// Returns approximate solutions by back-substitution through elimination ideals.
pub fn solve_elimination(generators: &[Polynomial], num_vars: usize) -> Vec<Vec<f64>> {
    if generators.is_empty() || num_vars == 0 {
        return vec![];
    }

    // Work through elimination ideals from last variable to first
    let gb = groebner_basis(generators, MonomialOrder::Lex);

    // Find polynomials in only the last variable
    let last_var_polys: Vec<&Polynomial> = gb
        .iter()
        .filter(|p| {
            p.terms.iter().all(|term| {
                term.exponents[..num_vars - 1].iter().all(|&e| e == 0)
            })
        })
        .collect();

    if last_var_polys.is_empty() {
        // No constraint on last variable, might be infinite solutions
        return vec![];
    }

    // Find approximate roots of univariate polynomial in last variable
    let last_poly = last_var_polys[0];
    let roots = find_approx_roots(last_poly, num_vars - 1, -10.0, 10.0, 200);

    // For each root, try to extend
    let mut solutions = Vec::new();
    for root_val in roots {
        let mut point = vec![0.0; num_vars];
        point[num_vars - 1] = root_val;

        // Back-substitute: for each variable from n-2 down to 0
        for var_idx in (0..num_vars - 1).rev() {
            let var_polys: Vec<&Polynomial> = gb
                .iter()
                .filter(|p| {
                    p.terms.iter().all(|term| {
                        term.exponents[..var_idx].iter().all(|&e| e == 0)
                            && term.exponents[var_idx..].iter().any(|&e| e > 0)
                    })
                })
                .collect();

            if let Some(vp) = var_polys.first() {
                // Try to find root for this variable given known values
                let candidates = find_approx_roots_partial(vp, var_idx, &point, -10.0, 10.0, 200);
                if let Some(&val) = candidates.first() {
                    point[var_idx] = val;
                }
            }
        }

        // Verify the solution
        let is_valid = generators.iter().all(|g| g.evaluate(&point).abs() < 1e-6);
        if is_valid {
            solutions.push(point);
        }
    }

    solutions
}

/// Find approximate roots of a univariate polynomial (in variable var_idx) by bisection.
fn find_approx_roots(poly: &Polynomial, var_idx: usize, lo: f64, hi: f64, steps: usize) -> Vec<f64> {
    let step = (hi - lo) / steps as f64;
    let mut roots = Vec::new();
    let mut prev_val = {
        let mut point = vec![0.0; poly.num_vars];
        point[var_idx] = lo;
        poly.evaluate(&point)
    };

    for i in 1..=steps {
        let x = lo + i as f64 * step;
        let mut point = vec![0.0; poly.num_vars];
        point[var_idx] = x;
        let val = poly.evaluate(&point);

        if val.abs() < 1e-10 {
            roots.push(x);
        } else if prev_val * val < 0.0 {
            // Sign change, bisect
            let root = bisect(poly, var_idx, x - step, x, 50);
            roots.push(root);
        }
        prev_val = val;
    }

    roots
}

/// Find approximate root given partial point values.
fn find_approx_roots_partial(
    poly: &Polynomial,
    var_idx: usize,
    point: &[f64],
    lo: f64,
    hi: f64,
    steps: usize,
) -> Vec<f64> {
    let step = (hi - lo) / steps as f64;
    let mut roots = Vec::new();
    let eval_at = |x: f64| -> f64 {
        let mut p = point.to_vec();
        p[var_idx] = x;
        poly.evaluate(&p)
    };

    let mut prev_val = eval_at(lo);

    for i in 1..=steps {
        let x = lo + i as f64 * step;
        let val = eval_at(x);

        if val.abs() < 1e-10 {
            roots.push(x);
        } else if prev_val * val < 0.0 {
            let mut a = x - step;
            let mut b = x;
            for _ in 0..50 {
                let mid = (a + b) / 2.0;
                let mid_val = eval_at(mid);
                if mid_val * eval_at(a) < 0.0 {
                    b = mid;
                } else {
                    a = mid;
                }
            }
            roots.push((a + b) / 2.0);
        }
        prev_val = val;
    }

    roots
}

fn bisect(poly: &Polynomial, var_idx: usize, lo: f64, hi: f64, iters: usize) -> f64 {
    let mut a = lo;
    let mut b = hi;
    for _ in 0..iters {
        let mid = (a + b) / 2.0;
        let mut point = vec![0.0; poly.num_vars];
        point[var_idx] = mid;
        let mid_val = poly.evaluate(&point);

        let mut point_a = vec![0.0; poly.num_vars];
        point_a[var_idx] = a;
        if mid_val * poly.evaluate(&point_a) < 0.0 {
            b = mid;
        } else {
            a = mid;
        }
    }
    (a + b) / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polynomial::PolynomialRing;

    fn ring2() -> PolynomialRing {
        PolynomialRing::new(2, MonomialOrder::Lex)
    }

    fn ring3() -> PolynomialRing {
        PolynomialRing::new(3, MonomialOrder::Lex)
    }

    #[test]
    fn test_elimination_ideal_basic() {
        let ring = ring2();
        // <x+y, x-y> with lex order: first elimination ideal should give y
        let f1 = ring.from_raw(vec![(1.0, vec![1, 0]), (1.0, vec![0, 1])]);
        let f2 = ring.from_raw(vec![(1.0, vec![1, 0]), (-1.0, vec![0, 1])]);
        let elim = elimination_ideal(&[f1, f2], 1);
        // Should give polynomial in y only
        assert!(!elim.is_empty());
        for p in &elim {
            assert!(p.terms.iter().all(|t| t.exponents[0] == 0));
        }
    }

    #[test]
    fn test_elimination_trivial() {
        let elim = elimination_ideal(&[], 1);
        assert!(elim.is_empty());
    }

    #[test]
    fn test_elimination_three_vars() {
        let ring = ring3();
        // <x + y + z, y + z>
        let f1 = ring.from_raw(vec![
            (1.0, vec![1, 0, 0]),
            (1.0, vec![0, 1, 0]),
            (1.0, vec![0, 0, 1]),
        ]);
        let f2 = ring.from_raw(vec![
            (1.0, vec![0, 1, 0]),
            (1.0, vec![0, 0, 1]),
        ]);
        let elim = elimination_ideal(&[f1, f2], 1);
        // After eliminating x, should have something in y, z
        for p in &elim {
            assert!(p.terms.iter().all(|t| t.exponents[0] == 0));
        }
    }

    #[test]
    fn test_solve_simple_system() {
        let ring = ring2();
        // x + y = 0, x - y = 0 => x = 0, y = 0
        let f1 = ring.from_raw(vec![(1.0, vec![1, 0]), (1.0, vec![0, 1])]);
        let f2 = ring.from_raw(vec![(1.0, vec![1, 0]), (-1.0, vec![0, 1])]);
        let sols = solve_elimination(&[f1, f2], 2);
        assert!(!sols.is_empty());
        for sol in &sols {
            assert!(sol[0].abs() < 0.1);
            assert!(sol[1].abs() < 0.1);
        }
    }

    #[test]
    fn test_can_extend() {
        let ring = ring2();
        let f1 = ring.from_raw(vec![(1.0, vec![1, 0]), (1.0, vec![0, 1])]); // x + y
        let f2 = ring.from_raw(vec![(1.0, vec![1, 0]), (-1.0, vec![0, 1])]); // x - y
        // x=0, y=0 can extend
        assert!(can_extend(&[f1, f2], 1, &[0.0]));
    }
}
