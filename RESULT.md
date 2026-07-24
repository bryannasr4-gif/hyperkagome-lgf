# RESULT — Hyperkagome lattice Green's function

**Addresses Varma–Monien, Phys. Rev. E 87, 032109 (2013)** ("we have currently found no way to
exactly solve the Eq. 16"; arXiv:1211.5666, submitted 2012).

> **Correction (July 2026).** The earlier headline — *not-a-symmetric-square ⇒ no elliptic closed form;
> SL₃-type Galois group* — is **retracted**. `Sym²(M)` has a rational solution, so `M`'s differential
> Galois group is **orthogonal** (`O(3,ℂ)`, `G° = SO(3,ℂ) ≅ PSL(2,ℂ)`) and `M` **is** projectively a
> symmetric square. Non-Liouvillianity survives. The Sym²/adjoint/intertwiner facts are certified in exact
> arithmetic (`numerics/certify_orthogonal.py`, `CERTIFICATE_orthogonal.txt`).
>
> **Resolution (July 2026, this repository).** The closed form is no longer merely *expected* — it is
> **proven and modular at level 30**. The second-order operator `V₂` (to which `M` is projectively
> `Sym²`-equivalent) is the **uniformizing ODE of the modular curve `X(Γ₀(30)⁺)`**: with the level-30
> eta quotient `u = [η(τ)η(6τ)η(10τ)η(15τ) / (η(2τ)η(3τ)η(5τ)η(30τ))]³`, the variable
> **`t = u/(u²+7u+1)`** (equivalently `1/t = u + 7 + 1/u`) **generates the genus-zero function field of
> `X(Γ₀(30)⁺)`**, and the
> Schwarzian identity `{τ,t} = 2·Q_V(t)` holds exactly. So the hyperkagome LGF is **modular at level
> `30 = 2·3·5`** — apparently the first lattice Green's function realised at a modular level with three
> distinct prime factors. Proven exactly (Ligozat + Atkin–Lehner + a pole-degree bound turning an 80-order
> series match into an identity) in `numerics/certify_modular.py` → `CERTIFICATE_modular.txt`; cross-vendor
> audited. This resolves the closed-form existence Varma & Monien conjectured. The weight-2 period
> `y₀ = Φ'/2` itself is exhibited in explicit closed form below — a depth-one quasimodular form on `Γ₀(30)⁺`
> twisted by the determinant character, **proven not** to be an eta quotient or (modular form)×(algebraic
> function) at any weight.

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
- **M is not a *literal* symmetric square** of a 2nd-order operator (exponent triples are not arithmetic
  progressions), certified at every singular point, with a runnable Sym² positive control returning the
  expected {0,1/2,1}. **But this does *not* exclude an elliptic form:** the AP test sees only
  function-multiplier and pullback equivalences, not operator homomorphisms.
- **Orthogonal differential Galois group.** `Sym²(M)` has an explicit rational solution
  `R(t) = −(1/272)(15t²+17t−8)² / (t²(t−1)²(4t−1)(5t−1)(9t−1))`; equivalently `M` is homomorphic to its
  adjoint by an order-2 intertwiner `T`. So the solution space carries a monodromy-invariant nondegenerate
  symmetric form ⇒ `G = O(3,ℂ)`, `G° = SO(3,ℂ) ≅ PSL(2,ℂ)`, and **M *is* projectively equivalent to the
  symmetric square of a 2nd-order operator** `V₂`; the det character of the monodromy = the quadratic
  character of the genus-one twist curve `v² = (1−4t)(1−5t)(1−9t)` (det = −1 at {1/9,1/5,1/4,∞}). Certified
  exactly in `numerics/certify_orthogonal.py` (Sym² membership margin 194 at t₀ = 1/2 and −1/3, made a proof
  by a Fuchs-relation budget with cap 109; exact intertwiner remainder identity; n = 2 Jordan; det character).
- **Modular parametrization — the closed form (proven).** `V₂` uniformizes `X(Γ₀(30)⁺)`:
  `t = u/(u²+7u+1)` with `u` the level-30 eta quotient above **generates the genus-zero function field of
  `X(Γ₀(30)⁺)`** (`1/t = u+7+1/u`; a *Weber-like function parametrization*, the term "Hauptmodul" being
  reserved for the `₂F₁` pullback `1728/j`), and `{τ,t} = 2·Q_V(t)` holds identically, where
  `Q_V = N(t)/[4t²(t−1)²(4t−1)²(5t−1)²(9t−1)²]`,
  `N = 24300t⁸−58860t⁷+73437t⁶−44294t⁵+15111t⁴−3160t³+407t²−30t+1`. The projective monodromy of `M` is the
  arithmetic lattice `Γ₀(30)⁺` (covolume `3π`, orbifold signature `(0;2,2,2,2,2;1 cusp)`, derived not
  assumed). The five order-2 points map to `t ∈ {1/9,1/5,1/4,1,∞}` (`t=1/9↔u=1`, `t=1/5↔u=−1`, the branch
  points of `u↦u+1/u`; `t=1/4,1,∞ ↔` roots of `u²+3u+1, u²+6u+1, u²+7u+1`) and the cusp to `t=0`. The twist
  curve `v²=(1−4t)(1−5t)(1−9t)` is `X(Γ₀(30)+2,3,6)`. ⇒ **hyperkagome LGF modular at level 30 = 2·3·5**.
  Proven exactly in `numerics/certify_modular.py` → `CERTIFICATE_modular.txt` (all checks PASS).
- **The weight-2 period `y₀ = Φ'/2` — explicit closed form.** `y₀ = [ρ₀(t)·W + ρ₁(t)·W'] / v` with
  `v = √((1−4t)(1−5t)(1−9t))`, `W = q·dt/dq`, `ρ₁ = (15t²+17t−8)/(30t(t−1))` and explicit rational `ρ₀`: a
  **weight-two, depth-one quasimodular** form on `Γ₀(30)⁺` twisted by the determinant character. The `W'` term
  is **provably essential** — `y₀` is not a (meromorphic modular form)×(algebraic function) at any weight or
  level, in particular **not an eta quotient** (proof: `G(V₂) = SL(2,ℂ)`, then a torus/unipotent argument in the
  Picard–Vessiot field). Certified exactly in `numerics/certify_y0.py` and `numerics/certify_y0_lemma.py`
  (→ `CERTIFICATE_y0.txt`, `CERTIFICATE_y0_lemma.txt`).
- **Non-Liouvillian: no algebraic or elementary closed form.** A genuine log at t=0 (repeated indicial
  exponent −1; exactly one log-free local solution of three, `numerics/certify_nonliouvillian.py`) excludes
  finite and imprimitive differential-Galois groups; equivalently `G° = SO(3,ℂ)` is simple hence
  non-solvable, so M is non-Liouvillian. This excludes only the algebraic/elementary cases — **not** the
  elliptic one (K itself is non-Liouvillian).
- **The Watson reduction (concrete realization; exactly verified).** `Φ(t) = 2/3 − (t/3) d/dt CT_θ log D₊`,
  `D₊ = (1−6t+3t²) − 2t²Σcos2θᵢ − 8t^{3/2}cosθ₁cosθ₂cosθ₃`: the hyperkagome LGF **is** the classical two-parameter
  generalized Watson integral (bcc lattice, 1st-neighbour weight `t^{3/2}` + 2nd-neighbour sc-shell weight `t²`) at
  spectral parameter `1−6t+3t²`. The `√t` cover = the `v²=(1−4t)(1−5t)(1−9t)` twist; `{1/9,1/5,1/4,1}` = band-map
  critical values (re-confirms the Riemann scheme). This anchors the level-30 modular parametrization above on the
  lattice side (both describe the same object). Verified in this repository
  (`numerics/verify_watson_reduction.py` → `CERTIFICATE_watson.txt`).
- **The E=1 reflection symmetry is structural:** hyperkagome is the line graph of a 3-regular
  *bipartite* premedial net Γ (8 vertices, 12 edges per cell); bipartite ⇒ adjacency spectrum symmetric
  about 0 ⇒ dispersive DOS symmetric about E=1 (line-graph shift d−2=1). Verified by a
  translation-consistent 2-colouring on finite blocks and by all odd central moments vanishing.

**Conditional-on-M note.** The structure theorem (irreducible; orthogonal Galois group; non-Liouvillian)
is a rigorous consequence *of M*, but M is established to the guess-and-verify standard (annihilation to
margin 55, independent rebuild), not by an unconditional creative-telescoping proof — which remains an
open item. (The orthogonality certificates themselves are exact, given M.) This is stated in the paper.

## Riemann scheme of M (Frobenius exponents)
- t=0            : {−1, −1, 0}   (repeated −1 ⇒ genuine log; all integers ⇒ single 3×3 Jordan block, n=2)
- t=1/9, 1/5, 1/4 : {−1/2, 0, 1}  (the −1/2 ⇒ √ van Hove DOS singularity at the band edge)
- t=1           : {−1, −1/2, 1/2} (principal van Hove point, E=0,2)
- t=∞           : {3/2, 2, 3} in the paper's convention (a solution of exponent ρ behaves as t^{−ρ});
                  equivalently {−3, −2, −3/2} as powers of t
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
- Not-*literal*-Sym² verdict from exact Frobenius exponents; SC positive control validates the test is
  non-vacuous.
- **Orthogonal structure (`numerics/certify_orthogonal.py`):** Sym²(M) rational solution R verified by
  exact series (194 relations at t₀ = 1/2 and −1/3, made a proof by a Fuchs-relation budget with cap 109;
  Gram matrix nondegenerate, rank 3, det ≠ 0);
  intertwiner remainder identity `rightremainder(M·T, adjoint(M)) = 0` exact over ℚ(t); indicial
  −64ρ(ρ+1)², log-free dim 1 ⇒ n = 2; Wronskian residues ⇒ det character −1 at {1/9,1/5,1/4,∞} ⇒ O(3,ℂ).
  All machinery validated on positive/negative controls (D³−4D self-dual, θ³ Jordan, etc.).

## Independent verification — status
- Lattice, moments (bit-identical to m₂₄₀), and operator independently reconstructed: AGREE.
- Irreducibility (over Q by complete exact enumeration; over Q̄ by Galois descent) and not-*literal*-Sym²:
  reproduced. The operator was reconstructed a second time independently, and the two reconstructions agree;
  the Sym²/adjoint homomorphism is certified here in exact arithmetic. An independent Maple `DFactor`
  / Magma factorization over Q̄ (operator provided in `numerics/M_maple.txt`, `M_magma.txt`, `M_sage.py`)
  and a `Homomorphisms`/`symmetric_power`+`ratsols` check are recommended external cross-checks.
- G(2) sign/exactness corrected: E=2 is singular; only Re G(1)=1/9 is exact.

Verdict: the structure theorem (irreducible order-3; not a *literal* symmetric square; Sym²(M) has a
rational solution ⇒ orthogonal Galois group G = O(3,ℂ), G° = SO(3,ℂ) ≅ PSL(2,ℂ) ⇒ projectively a symmetric
square; non-Liouvillian ⇒ no algebraic/elementary form) is sound and reproducible, conditional on M being the
minimal annihilator. The closed form is now **proven and modular**: `V₂` uniformizes `X(Γ₀(30)⁺)`, whose
function field is generated by `t = u/(u²+7u+1)` (`numerics/certify_modular.py`), so the hyperkagome LGF is **modular at level
30 = 2·3·5**, and the LGF is equivalently **the classical generalized bcc(1,2) Watson integral**
(`numerics/verify_watson_reduction.py`). The weight-2 period `y₀ = Φ'/2` is given in explicit depth-one
quasimodular closed form (`numerics/certify_y0.py`), **proven not** to reduce to an eta quotient or a
(modular form)×(algebraic function) at any weight (`numerics/certify_y0_lemma.py`). The one remaining open
item is an unconditional creative-telescoping proof that `M` is the minimal annihilator of the LGF.
