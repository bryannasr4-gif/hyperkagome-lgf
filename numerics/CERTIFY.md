# Optional independent cross-checks of the operator M

The repository's certification is **already complete in pure Python** (no SageMath required):

- `numerics/certify_factor.py` → `CERTIFICATE.txt` enumerates the complete finite set of Fuchs-admissible
  order-1 right factors of M and of adjoint(M) and tests each **exactly over ℚ** (0 found), computes the exact
  local exponents (including the degree-7 locus `{0,1,3}` by reduction mod p₇, and `t=∞` `{−3,−2,−3/2}`),
  proves a **genuine logarithmic solution at t=0** (exactly one log-free local solution of three), and hence
  concludes — by a Galois-descent argument — that **M is irreducible over ℚ̄(t)** and **not a symmetric square**.
- `numerics/certify_nonliouvillian.py` → `CERTIFICATE_nonliouvillian.txt` uses the same t=0 log to exclude finite
  and imprimitive Galois groups, i.e. **M is non-Liouvillian** (no algebraic/elementary closed form). The elliptic
  exclusion comes separately from *not-Sym²*.

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
# symmetric square: order 6 (not 3) confirms L is NOT a symmetric square:
S2 := symmetric_power(L, 2, [Dt, t]):  ratsols(S2, t);   # expect [] : no rational solution of Sym^2
formal_sol(L, [Dt, t], t = 0);         # expect a ln(t) term: the genuine logarithmic solution at t=0
```

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
- If any tool instead returns a lower-order factor, or if `ratsols`/`expsols` of the symmetric/exterior square are
  non-empty, that would contradict the pure-Python certificate and must be reconciled before any strong claim
  (this is exactly why the cross-check is offered).

Note: even with all cross-checks passing, the no-closed-form conclusion remains **conditional on M being the
minimal annihilator of the LGF** — established here to the guess-and-verify standard (annihilation margin 55), not
by an unconditional creative-telescoping proof. See `CT_SETUP.md` for that remaining step.
