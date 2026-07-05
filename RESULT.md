# RESULT — Hyperkagome lattice Green's function

**Addresses Varma–Monien, Phys. Rev. E 87, 032109 (2013)** ("we have currently found no way to
exactly solve the Eq. 16"; arXiv:1211.5666, submitted 2012).

## Main result
The local lattice Green's function G(z) of the 3D hyperkagome lattice (Na₄Ir₃O₈), after removing the
flat-band pole (1/3)/(z+2) and using the exact dispersive symmetry about E=1, has a generating
function Φ(t)=Σ ν_m t^m (t=(E−1)^{−2}) whose derivative Φ'(t) is annihilated by an **order-3
Picard–Fuchs operator M of degree 15** (the minimal operator L=M·d/dt has order 4).

- **M is irreducible over Q̄(t)** — certified by exact computation (`numerics/certify_factor.py`,
  `CERTIFICATE.txt`): a complete Fuchsian residue enumeration finds no order-1 right factor of M and
  none of adjoint(M) over Q (so no order-1 or order-2 right factor of M over Q), and a Galois-descent
  argument upgrades this to irreducibility over Q̄(t) using the genuine logarithmic solution at t=0
  (a hypothetical Q̄ right factor of order 1/2/3 is excluded respectively by the order-1 search, the
  adjoint/order-2 search, and the fact that an order-3 LCLM would force complete reducibility and hence
  no log). The enumeration uses the **exact** degree-7 exponents {0,1,3} (computed by reduction mod p₇),
  correcting an earlier floating-point read; infinity is confirmed regular singular ({−3,−2,−3/2}).
- **M is not a symmetric square** of a 2nd-order operator (exponent triples are not arithmetic
  progressions), certified at every singular point, with a runnable Sym² positive control returning the
  expected {0,1/2,1}. By the differential-Galois classification (Singer–Ulmer; van der Put–Singer),
  **⇒ no closed form in complete elliptic integrals** — in contrast to the SC/BCC/FCC lattices, which
  are Sym²(elliptic).
- **Non-Liouvillian (a separate exclusion): no algebraic or elementary closed form.** A genuine log
  at t=0 (repeated indicial exponent −1; exactly one log-free local solution of three,
  `numerics/certify_nonliouvillian.py`) excludes finite and imprimitive differential-Galois groups, so
  M is non-Liouvillian. NOTE this does **not** exclude the elliptic case (K itself is non-Liouvillian;
  the SC operator is irreducible + non-Liouvillian yet closed-form in K) — the elliptic exclusion comes
  solely from not-Sym². The two legs together: **no closed form in algebraic, elementary, OR elliptic
  functions.**
- **The E=1 reflection symmetry is structural:** hyperkagome is the line graph of a 3-regular
  *bipartite* premedial net Γ (8 vertices, 12 edges per cell); bipartite ⇒ adjacency spectrum symmetric
  about 0 ⇒ dispersive DOS symmetric about E=1 (line-graph shift d−2=1). Verified by a
  translation-consistent 2-colouring on finite blocks and by all odd central moments vanishing.

**Conditional-on-M note.** The no-closed-form conclusion is a rigorous consequence *of M*, but M is
established to the guess-and-verify standard (annihilation to margin 55, independent rebuild), not by an
unconditional creative-telescoping proof — which remains the open item. This is stated in the paper and
the emails.

## Riemann scheme of M (Frobenius exponents)
- t=0            : {−1, −1, 0}   (repeated −1 ⇒ the genuine logarithmic solution)
- t=1/9, 1/5, 1/4 : {−1/2, 0, 1}  (the −1/2 ⇒ √ van Hove DOS singularity at the band edge)
- t=1           : {−1, −1/2, 1/2} (principal van Hove point, E=0,2)
- t=∞           : {−3, −2, −3/2}
- p₇ locus (irreducible deg-7): {0, 1, 3}  (integer, apparent-type)

## Exact data (all verified)
- Flat band E=−2, weight exactly 1/3; band edges [−2,4]; dispersive DOS exactly symmetric about E=1.
- **Re G(1) = 1/9** — the headline exact special value: flat-band pole (1/3)/(1+2)=1/9 plus a dispersive
  principal value that vanishes by the E=1 symmetry. Confirmed numerically to ~1e-12 by symmetric BZ
  quadrature (`numerics/verify_specialvalues.py`); a naive Γ-inclusive grid gives a spurious few-percent
  offset, so use the off-Γ grid.
- **G(2): no clean exact value.** E=2 is the logarithmic van Hove point (t=1, a singular point of M),
  so Re G(2) is a singular boundary value (numerically ≈ +1/4, confirmable only to ~1–2%); it is *not*
  an exact special value. Re G(1)=1/9 is the sole headline exact value.
- Exact integer moments m₀..m₂₃₀ (saved); symmetry-reduced ν₀..ν₁₁₅.

## Certification summary
- Operator obtained by differential approximants (ore_algebra `guess`); **certified here in pure Python**.
- **Exact annihilation over Q of Φ'(t) through t^111: 112 relations vs 57 free coefficients, margin 55,
  0 residuals** (`numerics/strengthen_certification.py`, `numerics/verify.py`).
- Independent from-scratch rebuild of the lattice, moments (to m₂₄₀), and operator: bit-identical.
- Independent Bloch-Hamiltonian diagonalization reproduces the exact integer moments (a method disjoint
  from closed-walk enumeration).
- Not-Sym² verdict from exact Frobenius exponents; SC positive control validates the test is non-vacuous.

## Independent verification — status
- Lattice, moments (bit-identical to m₂₄₀), and operator independently reconstructed: AGREE.
- Irreducibility (over Q by complete exact enumeration; over Q̄ by Galois descent) and not-Sym²:
  reproduced. An independent Maple `DFactor` / Magma factorization over Q̄ (operator provided in
  `numerics/M_maple.txt`, `M_magma.txt`, `M_sage.py`) is a recommended external cross-check.
- G(2) sign/exactness corrected: E=2 is singular; only Re G(1)=1/9 is exact.

Verdict: headline result (irreducible order-3, not-Sym², non-Liouvillian ⇒ no algebraic/elementary/
elliptic closed form) is sound and reproducible, conditional on M being the minimal annihilator.
