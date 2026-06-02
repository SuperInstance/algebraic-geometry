# algebraic-geometry

Algebraic geometry in Rust. Varieties, ideals, and Groebner bases.

## What This Does

`algebraic-geometry` implements the core algorithms of computational algebraic geometry:

- **Multivariate polynomials** — arithmetic, evaluation, division, homogenization over ℚ (f64)
- **Monomial orderings** — lex, grlex, grevlex with full comparison
- **Polynomial ideals** — creation, membership testing, sum, product, intersection, equality
- **Gröbner bases** — Buchberger's algorithm with S-polynomials, reduction, and minimization
- **Affine and projective varieties** — point membership, dimension, Jacobians, singular points, union, intersection
- **Elimination theory** — elimination ideals, extension theorem, back-substitution solver
- **Bézout's theorem** — intersection counting and approximate intersection finding with Newton refinement
- **Polynomial constraint solving** — model constraints as polynomial systems and solve them

## Install

```bash
cargo add algebraic-geometry
```

Or in `Cargo.toml`:

```toml
[dependencies]
algebraic-geometry = "0.1"
```

Requires **Rust 2021 edition**. Dependencies: `serde`, `nalgebra`.

## Quick Start

```rust
use algebraic_geometry::*;
use algebraic_geometry::monomial_order::MonomialOrder;

let ring = PolynomialRing::new(2, MonomialOrder::Lex);

// x² - 1
let p = ring.from_raw(vec![
    (1.0, vec![2, 0]),   // x²
    (-1.0, vec![0, 0]),  // -1
]);

assert!(p.vanishes_at(&[1.0, 0.0], 1e-10));
```

## Testing

Run: `cargo test`

## License

MIT OR Apache-2.0
