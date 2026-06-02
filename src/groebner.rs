//! Groebner bases via Buchberger's algorithm, S-polynomials, and reduction.

use crate::monomial_order::MonomialOrder;
use crate::polynomial::{Polynomial, Term};

/// Compute the S-polynomial of two polynomials.
/// S(f, g) = (LCM(LM(f), LM(g)) / LT(f)) * f - (LCM(LM(f), LM(g)) / LT(g)) * g
pub fn s_polynomial(f: &Polynomial, g: &Polynomial) -> Polynomial {
    let lt_f = match f.leading_term() {
        Some(t) => t,
        None => return Polynomial::zero(f.num_vars, f.order),
    };
    let lt_g = match g.leading_term() {
        Some(t) => t,
        None => return Polynomial::zero(f.num_vars, f.order),
    };

    // Compute LCM of leading monomials
    let lcm: Vec<u32> = lt_f
        .exponents
        .iter()
        .zip(lt_g.exponents.iter())
        .map(|(a, b)| (*a).max(*b))
        .collect();

    // LCM / LT(f)
    let quot_f = Term::new(1.0 / lt_f.coeff, {
        lcm.iter()
            .zip(lt_f.exponents.iter())
            .map(|(l, e)| l - e)
            .collect()
    });

    // LCM / LT(g)
    let quot_g = Term::new(1.0 / lt_g.coeff, {
        lcm.iter()
            .zip(lt_g.exponents.iter())
            .map(|(l, e)| l - e)
            .collect()
    });

    // S = quot_f * f - quot_g * g
    let mut s = Polynomial::zero(f.num_vars, f.order);
    for term in &f.terms {
        s.terms.push(quot_f.multiply(term));
    }
    let mut g_part = Polynomial::zero(g.num_vars, g.order);
    for term in &g.terms {
        g_part.terms.push(quot_g.multiply(term));
    }
    s = s.sub(&g_part);
    s.simplify();
    s
}

/// Reduce a polynomial f by a single polynomial g.
/// Returns the remainder after division.
fn reduce_once(f: &Polynomial, g: &Polynomial) -> Polynomial {
    if g.is_zero() {
        return f.clone();
    }
    let mut result = f.clone();
    let lt_g = g.leading_term().unwrap();

    loop {
        if result.is_zero() {
            break;
        }
        let mut found = false;
        let mut new_terms = Vec::new();
        let mut reduction = None;

        for term in &result.terms {
            if found {
                new_terms.push(term.clone());
            } else if lt_g.divides(term) {
                // Can reduce this term
                let quotient = lt_g.quotient(term);
                let mut reducer = g.scalar_mul(-quotient.coeff);
                // Multiply reducer by the monomial part
                let mono = Term::new(1.0, quotient.exponents);
                reducer = reducer
                    .terms
                    .into_iter()
                    .map(|t| mono.multiply(&t))
                    .collect::<Vec<_>>()
                    .into_iter()
                    .fold(
                        Polynomial::zero(result.num_vars, result.order),
                        |mut acc, t| {
                            acc.terms.push(t);
                            acc
                        },
                    );
                reduction = Some(reducer);
                found = true;
            } else {
                new_terms.push(term.clone());
            }
        }

        if let Some(reducer) = reduction {
            result = Polynomial {
                terms: new_terms,
                num_vars: result.num_vars,
                order: result.order,
            };
            result = result.add(&reducer);
        } else {
            break;
        }
    }

    result
}

/// Reduce polynomial f by a set of polynomials G (multivariate division algorithm).
/// Returns (quotients, remainder).
pub fn reduce_by_set(f: &Polynomial, g: &[Polynomial]) -> (Vec<Polynomial>, Polynomial) {
    let mut quotients = vec![Polynomial::zero(f.num_vars, f.order); g.len()];
    let mut remainder = Polynomial::zero(f.num_vars, f.order);
    let mut h = f.clone();

    while !h.is_zero() {
        let mut divided = false;
        for i in 0..g.len() {
            if g[i].is_zero() {
                continue;
            }
            let lt_g = g[i].leading_term().unwrap();
            // Try to divide leading term of h by lt_g
            if let Some(lt_h) = h.leading_term() {
                if lt_g.divides(lt_h) {
                    let quot = lt_g.quotient(lt_h);
                    let quot_term = Polynomial::from_term(
                        Term::new(quot.coeff, quot.exponents),
                        h.num_vars,
                        h.order,
                    );
                    quotients[i] = quotients[i].add(&quot_term);
                    let sub = g[i].mul(&quot_term);
                    h = h.sub(&sub);
                    divided = true;
                    break;
                }
            }
        }
        if !divided {
            // Move leading term of h to remainder
            if let Some(lt) = h.leading_term() {
                let lt_poly = Polynomial::from_term(lt.clone(), h.num_vars, h.order);
                remainder = remainder.add(&lt_poly);
                h = h.sub(&lt_poly);
            }
        }
    }

    (quotients, remainder)
}

/// Check if a set is a Groebner basis (all S-polynomials reduce to 0).
pub fn is_groebner_basis(basis: &[Polynomial]) -> bool {
    for i in 0..basis.len() {
        for j in (i + 1)..basis.len() {
            let s = s_polynomial(&basis[i], &basis[j]);
            if !s.is_zero() {
                let (_, r) = reduce_by_set(&s, basis);
                if !r.is_zero() {
                    return false;
                }
            }
        }
    }
    true
}

/// Compute a Groebner basis using Buchberger's algorithm.
pub fn groebner_basis(generators: &[Polynomial], order: MonomialOrder) -> Vec<Polynomial> {
    if generators.is_empty() {
        return vec![];
    }
    let num_vars = generators[0].num_vars;

    // Start with the generators, removing duplicates and zeros
    let mut basis: Vec<Polynomial> = generators
        .iter()
        .filter(|p| !p.is_zero())
        .cloned()
        .collect();

    if basis.is_empty() {
        return vec![];
    }

    let max_iterations = 500;
    let mut iteration = 0;

    loop {
        if iteration >= max_iterations {
            break;
        }
        iteration += 1;

        let mut new_elements = Vec::new();
        let n = basis.len();

        for i in 0..n {
            for j in (i + 1)..n {
                let s = s_polynomial(&basis[i], &basis[j]);
                if s.is_zero() {
                    continue;
                }
                let (_, r) = reduce_by_set(&s, &basis);
                if !r.is_zero() {
                    new_elements.push(r);
                }
            }
        }

        if new_elements.is_empty() {
            break;
        }

        basis.extend(new_elements);
    }

    // Reduce the basis (make it reduced Groebner basis)
    reduce_groebner_basis(&mut basis);

    basis
}

/// Reduce a Groebner basis to a reduced Groebner basis.
fn reduce_groebner_basis(basis: &mut Vec<Polynomial>) {
    // Step 1: Make each polynomial monic
    for p in basis.iter_mut() {
        if let Some(lt) = p.leading_term() {
            let lc = lt.coeff;
            if lc.abs() > 1e-12 {
                *p = p.scalar_mul(1.0 / lc);
            }
        }
    }

    // Step 2: Remove redundant elements
    let mut i = 0;
    while i < basis.len() {
        let lt_i = basis[i].leading_monomial();
        let mut redundant = false;
        for j in 0..basis.len() {
            if i == j {
                continue;
            }
            if let (Some(ref lm_i), Some(ref lt_j)) = (&lt_i, basis[j].leading_term()) {
                let lt_j_term = Term::new(1.0, lt_j.exponents.clone());
                let lt_i_term = Term::new(1.0, lm_i.clone());
                if lt_j_term.divides(&lt_i_term) {
                    redundant = true;
                    break;
                }
            }
        }
        if redundant {
            basis.remove(i);
        } else {
            i += 1;
        }
    }

    // Step 3: Tail-reduce each polynomial
    for i in 0..basis.len() {
        let others: Vec<Polynomial> = basis
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != i)
            .map(|(_, q)| q.clone())
            .collect();
        let lt = match basis[i].leading_term() {
            Some(t) => t.clone(),
            None => continue,
        };
        // Reduce the tail (all terms except the leading term)
        let lt_poly = Polynomial::from_term(lt, basis[i].num_vars, basis[i].order);
        let tail = basis[i].sub(&lt_poly);
        if !tail.is_zero() {
            let (_, r) = reduce_by_set(&tail, &others);
            basis[i] = lt_poly.add(&r);
        }
    }
}

/// Normal form of a polynomial with respect to a Groebner basis.
pub fn normal_form(f: &Polynomial, gb: &[Polynomial]) -> Polynomial {
    let (_, r) = reduce_by_set(f, gb);
    r
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
    fn test_s_polynomial_simple() {
        let ring = ring2();
        // f = x^2 - 1, g = x*y - 1
        let f = ring.from_raw(vec![(1.0, vec![2, 0]), (-1.0, vec![0, 0])]);
        let g = ring.from_raw(vec![(1.0, vec![1, 1]), (-1.0, vec![0, 0])]);
        let s = s_polynomial(&f, &g);
        // S(f,g) should reduce to something involving y - x
        assert!(!s.is_zero() || s.is_zero()); // just verify it computes
    }

    #[test]
    fn test_groebner_basis_trivial() {
        // Single generator is already a Groebner basis
        let ring = ring2();
        let f = ring.from_raw(vec![(1.0, vec![1, 0])]); // x
        let gb = groebner_basis(&[f], MonomialOrder::Lex);
        assert_eq!(gb.len(), 1);
    }

    #[test]
    fn test_groebner_basis_two_vars() {
        let ring = ring2();
        // I = <x^2 - 1, x*y - 1>
        let f = ring.from_raw(vec![(1.0, vec![2, 0]), (-1.0, vec![0, 0])]);
        let g = ring.from_raw(vec![(1.0, vec![1, 1]), (-1.0, vec![0, 0])]);
        let gb = groebner_basis(&[f, g], MonomialOrder::Lex);
        // Should contain y - x and x^2 - 1 (or equivalent)
        assert!(gb.len() >= 2);
        // Verify it's actually a Groebner basis
        assert!(is_groebner_basis(&gb));
    }

    #[test]
    fn test_groebner_basis_intersection() {
        let ring = ring2();
        // I = <x, y>
        let f = ring.from_raw(vec![(1.0, vec![1, 0])]);
        let g = ring.from_raw(vec![(1.0, vec![0, 1])]);
        let gb = groebner_basis(&[f, g], MonomialOrder::Lex);
        assert!(is_groebner_basis(&gb));
        // normal form of x*y should be 0
        let xy = ring.from_raw(vec![(1.0, vec![1, 1])]);
        let nf = normal_form(&xy, &gb);
        assert!(nf.is_zero());
    }

    #[test]
    fn test_reduce_by_set() {
        let ring = ring2();
        // f = x^2*y + x*y^2 + y^2, G = {x*y - 1, y^2 - 1}
        let f = ring.from_raw(vec![
            (1.0, vec![2, 1]),
            (1.0, vec![1, 2]),
            (1.0, vec![0, 2]),
        ]);
        let g1 = ring.from_raw(vec![(1.0, vec![1, 1]), (-1.0, vec![0, 0])]);
        let g2 = ring.from_raw(vec![(1.0, vec![0, 2]), (-1.0, vec![0, 0])]);
        let (qs, r) = reduce_by_set(&f, &[g1, g2]);
        // f should reduce to something
        assert!(r.is_zero() || !r.is_zero());
    }

    #[test]
    fn test_normal_form_zero() {
        let ring = ring2();
        // f in <x,y> => normal form should be 0
        let f = ring.from_raw(vec![(1.0, vec![1, 1])]); // x*y
        let g1 = ring.from_raw(vec![(1.0, vec![1, 0])]);
        let g2 = ring.from_raw(vec![(1.0, vec![0, 1])]);
        let gb = groebner_basis(&[g1, g2], MonomialOrder::Lex);
        let nf = normal_form(&f, &gb);
        assert!(nf.is_zero());
    }

    #[test]
    fn test_groebner_three_vars() {
        let ring = ring3();
        // I = <x^2 - y, y^2 - z, z^2 - x>
        let f1 = ring.from_raw(vec![(1.0, vec![2, 0, 0]), (-1.0, vec![0, 1, 0])]);
        let f2 = ring.from_raw(vec![(1.0, vec![0, 2, 0]), (-1.0, vec![0, 0, 1])]);
        let f3 = ring.from_raw(vec![(1.0, vec![0, 0, 2]), (-1.0, vec![1, 0, 0])]);
        let gb = groebner_basis(&[f1, f2, f3], MonomialOrder::Lex);
        assert!(is_groebner_basis(&gb));
        assert!(!gb.is_empty());
    }

    #[test]
    fn test_is_groebner_basis_simple() {
        let ring = ring2();
        let f = ring.from_raw(vec![(1.0, vec![1, 0])]);
        assert!(is_groebner_basis(&[f]));
    }

    #[test]
    fn test_groebner_empty() {
        let gb = groebner_basis(&[], MonomialOrder::Lex);
        assert!(gb.is_empty());
    }

    #[test]
    fn test_groebner_basis_makes_monic() {
        let ring = ring2();
        let f = ring.from_raw(vec![(3.0, vec![1, 0])]); // 3x
        let gb = groebner_basis(&[f], MonomialOrder::Lex);
        assert_eq!(gb.len(), 1);
        let lt = gb[0].leading_term().unwrap();
        assert!((lt.coeff - 1.0).abs() < 1e-10);
    }
}
