# Optional independent cross-checks of the operator M

The repository's certification is **already complete in pure Python** (no SageMath required):

- `numerics/certify_factor.py` → `CERTIFICATE.txt` enumerates the complete finite set of Fuchs-admissible
  order-1 right factors of M and of adjoint(M) and tests each **exactly over ℚ** (0 found), computes the exact
  local exponents (including the degree-7 locus `{0,1,3}` by reduction mod p₇, and `t=∞`, printed as the
  *powers of t* `{−3,−2,−3/2}`, i.e. `y ~ t^−3, t^−2, t^−3/2`; in the Riemann-scheme convention of the paper,
  where an exponent ρ at ∞ means `y ~ t^−ρ`, the same three solutions are written `{3/2,2,3}` — that is the
  convention under which the Fuchs relation `Σ exponents = n(n−1)/2·(s−2) = 33` holds),
  proves a **genuine logarithmic solution at t=0** (exactly one log-free local solution of three), and hence
  concludes — by a Galois-descent argument — that **M is irreducible over ℚ̄(t)** and **not a *literal* symmetric
  square**.
- `numerics/certify_orthogonal.py` → `CERTIFICATE_orthogonal.txt` shows that **Sym²(M) has a rational solution**
  and that **M is homomorphic to its adjoint** (order-2 intertwiner), so the differential Galois group is
  **orthogonal**: `G = O(3,C)`, `G⁰ = SO(3,C) ≅ PSL(2,C)`, and **M IS projectively a symmetric square** ⇒ an
  elliptic/modular closed form is *expected*. (These facts are certified here in exact arithmetic.)
- `numerics/certify_nonliouvillian.py` → `CERTIFICATE_nonliouvillian.txt` uses the same t=0 log to exclude finite
  and imprimitive Galois groups, i.e. **M is non-Liouvillian** (no algebraic/elementary closed form). This does
  **not** exclude an elliptic form — the earlier "not-Sym² ⇒ no elliptic" claim is **retracted**.

None of the below is required. They are **independent cross-checks in other computer-algebra systems**, useful if
a referee wants confirmation from a certified decision procedure. All use the **same operator M**, provided in
`M_maple.txt`, `M_magma.txt`, and `M_sage.py` (each defines `c0..c3` and `L = c3·D³ + c2·D² + c1·D + c0`).

## Route A — Maple (`DFactor` handles order ≤ 4 well)

```maple
with(DEtools):
# paste numerics/M_maple.txt (defines c0,c1,c2,c3 and L):
F := DFactor(L, [Dt, t]);         # a single factor (= L) => IRREDUCIBLE over Q(t)
nops([F]);                        # 1  => irreducible ; >1 => reducible, print the factors
eigenring(L, [Dt, t]);            # cross-check: eigenring = scalars (dim 1) <=> irreducible
# L is not a LITERAL symmetric square, but it is homomorphic to its adjoint (orthogonal Galois group):
Homomorphisms(adjoint(L), L, [Dt, t]);   # expect a NONZERO order-2 intertwiner T
S2 := symmetric_power(L, 2, [Dt, t]):  ratsols(S2, t);
    # expect a NONZERO rational solution, R = (15t^2+17t-8)^2/(t^2(t-1)^2(4t-1)(5t-1)(9t-1))
    # up to the arbitrary overall constant ratsols happens to normalize to
    # (=> Galois group in O(3,C); confirms certify_orthogonal.py)
formal_sol(L, [Dt, t], t = 0);         # expect a ln(t) term: the genuine logarithmic solution at t=0
```

### The intertwiner with `Sym^2(V2)`: conjugate by the twist first

Searching for it directly returns **nothing**, and that is forced, not a Maple limitation: the determinant
character of `L` is the quadratic character of `v^2 = (1-4t)(1-5t)(1-9t)`, while `Sym^2(V2)` in projective
normal form is unimodular, so the two are non-isomorphic (`numerics/certify_intertwiner.py`). Conjugating by
`v` removes exactly that obstruction, and `v'/v` is rational so the conjugate stays in `Q(t)<Dt>`:

```maple
QV := (1/4)*(24300*t^8 - 58860*t^7 + 73437*t^6 - 44294*t^5 + 15111*t^4 - 3160*t^3
      + 407*t^2 - 30*t + 1)/t^2/(t-1)^2/(4*t-1)^2/(5*t-1)^2/(9*t-1)^2:
V2 := Dt^2 + QV:
v  := (1-4*t)^(1/2)*(1-5*t)^(1/2)*(1-9*t)^(1/2):
Homomorphisms(symmetric_power(V2, 2), L, [Dt, t]);              # [] -- empty, necessarily
Homomorphisms(symmetric_power(V2, 2),
              mult(v, L, 1/v, [Dt, t]), [Dt, t]);               # the ORDER-ONE T = rho0 + rho1*Dt
```

with `rho1 = (15t^2+17t-8)/(30 t (t-1))` and `rho0` as in Eq. (12) of the paper. The obstruction is a character
of order two, so it cannot survive an even tensor construction: `det Sym^2 = (det)^4`, which is why
`Homomorphisms(symmetric_power(V2,4), symmetric_power(L,2))` **is** nonempty without any conjugation.
(Recipe due to J.-M. Maillard, private communication, July-August 2026.)

Read-off: `DFactor` returns one operator ⇒ irreducible over ℚ(t) (and, with the exact enumeration over ℚ̄ plus
the t=0 log already in `certify_factor.py`, over ℚ̄(t)). A `ln(t)` term at t=0 confirms the logarithm.

## Route B — Magma (`Factorisation`, a certified decision procedure)

```magma
// paste numerics/M_magma.txt (defines F,R,D, c0..c3, L):
fac := Factorisation(L);
#fac;                          // 1 => irreducible order-3 operator
[Order(f[1]) : f in fac];
```

## Route C — Sage + ore_algebra (free; use a fresh conda env)

```bash
conda create -n sage2 -c conda-forge sage python=3.11 -y
conda activate sage2 && sage -pip install ore_algebra
sage numerics/M_sage.py        # builds the same L and calls L.factor(); one order-3 factor => irreducible
```

## Interpretation

- A single order-3 factor from `DFactor` / `Factorisation` / `factor()` **independently confirms** the
  irreducibility already certified in pure Python.
- `Homomorphisms(adjoint(L), L)` should return a **nonzero** order-2 intertwiner, and `ratsols(symmetric_power(L,2))`
  a **nonzero** rational solution — both **expected** (they establish the orthogonal Galois group). If instead they
  came back empty, *that* would contradict the pure-Python `certify_orthogonal.py` certificate and must be
  reconciled. If any factorizer returns a lower-order factor, that would contradict the irreducibility certificate.

Note: even with all cross-checks passing, the structure theorem remains **conditional on M being the minimal
annihilator of the LGF** — established here to the guess-and-verify standard (annihilation margin 55), not by an
unconditional creative-telescoping proof. See `CT_SETUP.md` for that remaining step.
