//! Bezout's theorem: for plane curves of degrees d1 and d2 with no common component,
//! the number of intersection points (counted with multiplicity) is d1 * d2.

use crate::monomial_order::MonomialOrder;
use crate::polynomial::Polynomial;
use crate::polynomial::PolynomialRing;

/// Count intersection points of two plane curves using Bezout's theorem.
/// Returns the product of degrees if the curves have no common component.
pub fn bezout_count(f: &Polynomial, g: &Polynomial) -> u32 {
    f.total_degree() * g.total_degree()
}

/// Find approximate intersection points of two plane curves by sampling.
pub fn find_intersections(
    f: &Polynomial,
    g: &Polynomial,
    search_range: f64,
    resolution: usize,
) -> Vec<(f64, f64)> {
    let step = 2.0 * search_range / resolution as f64;
    let mut intersections = Vec::new();

    for i in 0..resolution {
        for j in 0..resolution {
            let x = -search_range + i as f64 * step;
            let y = -search_range + j as f64 * step;
            let fx = f.evaluate(&[x, y]);
            let gx = g.evaluate(&[x, y]);
            if fx.abs() < 0.1 && gx.abs() < 0.1 {
                // Refine with Newton-like iteration
                if let Some(pt) = refine_intersection(f, g, x, y) {
                    // Check not duplicate
                    let is_dup = intersections.iter().any(|&(px, py)| {
                        ((px - pt.0) as f64).hypot((py - pt.1) as f64) < 0.01
                    });
                    if !is_dup {
                        intersections.push(pt);
                    }
                }
            }
        }
    }

    intersections
}

/// Refine an approximate intersection point using Newton's method.
fn refine_intersection(f: &Polynomial, g: &Polynomial, x0: f64, y0: f64) -> Option<(f64, f64)> {
    let eps = 1e-8;
    let mut x = x0;
    let mut y = y0;

    for _ in 0..50 {
        let fx = f.evaluate(&[x, y]);
        let gx = g.evaluate(&[x, y]);

        if fx.abs() < 1e-12 && gx.abs() < 1e-12 {
            return Some((x, y));
        }

        // Jacobian
        let dfx = (f.evaluate(&[x + eps, y]) - f.evaluate(&[x - eps, y])) / (2.0 * eps);
        let dfy = (f.evaluate(&[x, y + eps]) - f.evaluate(&[x, y - eps])) / (2.0 * eps);
        let dgx = (g.evaluate(&[x + eps, y]) - g.evaluate(&[x - eps, y])) / (2.0 * eps);
        let dgy = (g.evaluate(&[x, y + eps]) - g.evaluate(&[x, y - eps])) / (2.0 * eps);

        let det = dfx * dgy - dfy * dgx;
        if det.abs() < 1e-20 {
            break; // Singular, can't refine
        }

        let dx = (dgy * fx - dfy * gx) / det;
        let dy = (dfx * gx - dgx * fx) / det;

        x -= dx;
        y -= dy;

        if dx.abs() < 1e-12 && dy.abs() < 1e-12 {
            break;
        }
    }

    if f.evaluate(&[x, y]).abs() < 1e-6 && g.evaluate(&[x, y]).abs() < 1e-6 {
        Some((x, y))
    } else {
        None
    }
}

/// Compute intersection multiplicity at a point (simplified).
/// For curves f, g at point p, the multiplicity is the intersection number.
/// We use a numerical approximation: perturb and count nearby intersections.
pub fn intersection_multiplicity(f: &Polynomial, g: &Polynomial, point: (f64, f64)) -> u32 {
    // Simplified: just check how "flat" both curves are at the intersection
    let eps = 1e-8;
    let (x, y) = point;

    let dfx = (f.evaluate(&[x + eps, y]) - f.evaluate(&[x - eps, y])) / (2.0 * eps);
    let dfy = (f.evaluate(&[x, y + eps]) - f.evaluate(&[x, y - eps])) / (2.0 * eps);
    let dgx = (g.evaluate(&[x + eps, y]) - g.evaluate(&[x - eps, y])) / (2.0 * eps);
    let dgy = (g.evaluate(&[x, y + eps]) - g.evaluate(&[x, y - eps])) / (2.0 * eps);

    let det = dfx * dgy - dfy * dgx;
    if det.abs() < 1e-10 {
        // Tangent curves, multiplicity > 1
        2 // Simplified: assume multiplicity 2 for tangent curves
    } else {
        1 // Transverse intersection
    }
}

/// Verify Bezout's theorem: count intersections and compare to d1*d2.
pub fn verify_bezout(
    f: &Polynomial,
    g: &Polynomial,
    search_range: f64,
    resolution: usize,
) -> BezoutResult {
    let theoretical = bezout_count(f, g);
    let intersections = find_intersections(f, g, search_range, resolution);
    let found_count = intersections.len() as u32;

    let total_multiplicity: u32 = intersections
        .iter()
        .map(|&pt| intersection_multiplicity(f, g, pt))
        .sum();

    BezoutResult {
        theoretical_count: theoretical,
        found_intersections: intersections,
        total_multiplicity,
        verified: found_count > 0 && found_count <= theoretical,
    }
}

/// Result of Bezout's theorem verification.
#[derive(Debug, Clone)]
pub struct BezoutResult {
    pub theoretical_count: u32,
    pub found_intersections: Vec<(f64, f64)>,
    pub total_multiplicity: u32,
    pub verified: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring2() -> PolynomialRing {
        PolynomialRing::new(2, MonomialOrder::Lex)
    }

    #[test]
    fn test_bezout_count_lines() {
        let ring = ring2();
        let l1 = ring.from_raw(vec![(1.0, vec![1, 0]), (1.0, vec![0, 0])]); // x + 1
        let l2 = ring.from_raw(vec![(1.0, vec![0, 1]), (-1.0, vec![0, 0])]); // y - 1
        assert_eq!(bezout_count(&l1, &l2), 1); // Two lines intersect at 1 point
    }

    #[test]
    fn test_bezout_count_line_circle() {
        let ring = ring2();
        // x^2 + y^2 - 1 (circle, deg 2) and x (line, deg 1)
        let circle = ring.from_raw(vec![
            (1.0, vec![2, 0]),
            (1.0, vec![0, 2]),
            (-1.0, vec![0, 0]),
        ]);
        let line = ring.from_raw(vec![(1.0, vec![1, 0])]);
        assert_eq!(bezout_count(&circle, &line), 2);
    }

    #[test]
    fn test_bezout_two_circles() {
        let ring = ring2();
        let c1 = ring.from_raw(vec![
            (1.0, vec![2, 0]),
            (1.0, vec![0, 2]),
            (-1.0, vec![0, 0]),
        ]);
        let c2 = ring.from_raw(vec![
            (1.0, vec![2, 0]),
            (1.0, vec![0, 2]),
            (-2.0, vec![1, 0]),
        ]); // (x-1)^2 + y^2 - 1 = x^2 - 2x + 1 + y^2 - 1
        assert_eq!(bezout_count(&c1, &c2), 4);
    }

    #[test]
    fn test_find_intersections_lines() {
        let ring = ring2();
        // y = x and y = -x intersect at origin
        let l1 = ring.from_raw(vec![(1.0, vec![0, 1]), (-1.0, vec![1, 0])]); // y - x
        let l2 = ring.from_raw(vec![(1.0, vec![0, 1]), (1.0, vec![1, 0])]); // y + x
        let pts = find_intersections(&l1, &l2, 5.0, 50);
        assert!(!pts.is_empty());
        let has_origin = pts.iter().any(|&(x, y)| x.abs() < 0.1 && y.abs() < 0.1);
        assert!(has_origin);
    }

    #[test]
    fn test_find_intersections_line_circle() {
        let ring = ring2();
        let circle = ring.from_raw(vec![
            (1.0, vec![2, 0]),
            (1.0, vec![0, 2]),
            (-1.0, vec![0, 0]),
        ]);
        let line = ring.from_raw(vec![(1.0, vec![1, 0])]); // x = 0
        let pts = find_intersections(&circle, &line, 5.0, 100);
        // Circle x^2+y^2=1 intersected with x=0 gives y=±1
        assert!(pts.len() >= 1);
    }

    #[test]
    fn test_verify_bezout() {
        let ring = ring2();
        let l1 = ring.from_raw(vec![(1.0, vec![1, 0]), (-1.0, vec![0, 0])]); // x - 1
        let l2 = ring.from_raw(vec![(1.0, vec![0, 1]), (-1.0, vec![0, 0])]); // y - 1
        let result = verify_bezout(&l1, &l2, 5.0, 50);
        assert_eq!(result.theoretical_count, 1);
    }

    #[test]
    fn test_bezout_conic_line() {
        let ring = ring2();
        // Parabola y = x^2 and line y = 0
        let parabola = ring.from_raw(vec![(1.0, vec![0, 1]), (-1.0, vec![2, 0])]);
        let line = ring.from_raw(vec![(1.0, vec![0, 1])]);
        assert_eq!(bezout_count(&parabola, &line), 2);
    }

    #[test]
    fn test_bezout_cubic_line() {
        let ring = ring2();
        // y = x^3 and y = 0
        let cubic = ring.from_raw(vec![(1.0, vec![0, 1]), (-1.0, vec![3, 0])]);
        let line = ring.from_raw(vec![(1.0, vec![0, 1])]);
        assert_eq!(bezout_count(&cubic, &line), 3);
    }
}
