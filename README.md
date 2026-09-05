# Hyperkagome lattice Green's function: modular uniformization at level 30

This repository accompanies the manuscript *"Lattice Green's function of the hyperkagome lattice:
modular uniformization at level 30 from an orthogonal differential Galois group"*, by **Bryan Nasr
and Jean-Marie Maillard** ([`paper/main.pdf`](paper/main.pdf)). It contains the exact operator, the
data it was built from, and self-contained scripts that reproduce every certification.

**Preprint:** [arXiv:2608.28141](https://arxiv.org/abs/2608.28141)
(`math-ph`, cross-listed `cond-mat.stat-mech` and `math.NT`),
[doi:10.48550/arXiv.2608.28141](https://doi.org/10.48550/arXiv.2608.28141). The code in this
repository is MIT-licensed; the preprint is distributed under the arXiv.org perpetual
non-exclusive licence.

> **Correction notice (July 2026).** An earlier version of this work claimed that `M`, being irreducible
> and not a symmetric square, has *no closed form in complete elliptic integrals* and an "SL₃-type"
> differential Galois group. **That claim was wrong and has been retracted.** `Sym²(M)` has a rational
> solution, so `M`'s Galois group is orthogonal (`O(3,ℂ)`, identity component `SO(3,ℂ) ≅ PSL(2,ℂ)`) and
> `M` **is** projectively equivalent to the symmetric square of a second-order operator `V₂`. The
> homomorphism-to-adjoint, the rational solution of `Sym²(M)`, and the explicit intertwiner establish the
> orthogonal structure; they are certified in exact arithmetic
> (`numerics/certify_orthogonal.py`). Non-Liouvillianity (no algebraic/elementary form) still holds.
>
> **Resolution (July 2026).** The closed form is no longer merely *expected* — it is **proven and modular
> at level 30**. `V₂` is the **uniformizing ODE of the modular curve `X(Γ₀(30)⁺)`**: with the level-30
> eta quotient `u = [η(τ)η(6τ)η(10τ)η(15τ) / (η(2τ)η(3τ)η(5τ)η(30τ))]³`, the variable **`t = u/(u²+7u+1)`**
> (i.e. `1/t = u+7+1/u`) **generates the genus-zero function field of `X(Γ₀(30)⁺)`**, and `{τ,t} = 2·Q_V(t)` holds exactly. So the
> hyperkagome LGF is **modular at level `30 = 2·3·5`** — apparently the first lattice Green's function at a
> modular level with three distinct prime factors. Proven exactly in
> [`numerics/certify_modular.py`](numerics/certify_modular.py) → `CERTIFICATE_modular.txt`.
>
> **Prior art, and what is new (August 2026).** That uniformizing equation is itself not new: it is the
> row `Γ₀(30)⁺` of the table of genus-zero groups in **B. H. Lian and S.-T. Yau**, *Mirror maps,
> modular relations and hypergeometric series I* (arXiv:hep-th/9507151), where it arises as the
> Picard–Fuchs operator of a degenerating family of algebraic K3 surfaces, and its Schwarzian potential
> is the tabulated Conway–Norton class `30B` of Lian–Wiczer (arXiv:math/0611291). **What is new here is
> the identification of a lattice Green's function with that operator, together with the proof**: the
> Schwarzian identity is established by an a priori pole-degree bound rather than by a series match, and
> therefore also proves the tabulated entries, which are exhibited rather than bounded. Both
> identifications are verified in exact rational arithmetic in
> [`numerics/certify_tabulated.py`](numerics/certify_tabulated.py) → `CERTIFICATE_tabulated.txt`. The
> comparison needs the quadratic-differential Jacobian: a Schwarzian `Q` is not a function, and a naive
> substitution differs from the tabulated entry by exactly `(1+4z)⁴`.

## Result

Varma & Monien (*Lattice Green's functions for kagome, diced and hyperkagome lattices*,
[Phys. Rev. E **87**, 032109 (2013)](https://doi.org/10.1103/PhysRevE.87.032109); arXiv:1211.5666)
reduced the 3D hyperkagome (Na₄Ir₃O₈) density of states to a threefold integral (their Eq. 16) and
wrote that they "found no way to exactly solve" it.

Working from exact lattice moments, we determine that the hyperkagome Green's function — after
removing the flat-band pole and using the exact reflection symmetry about `E = 1` — is annihilated
by an **irreducible order-3 Picard–Fuchs operator `M` of degree 15**, and that `M` is **not a
*literal* symmetric square** of a second-order operator (its Frobenius exponent triples are not the
arithmetic progressions a symmetric square would force). The order-4 operator guessed from the
moments factors as `L₄ = M·d/dt`, and that is its **unique** factorization: `L₄` is *not*
`LCLM(N, d/dt)` for any order-3 `N`, so the constant solution is not a direct summand and `M` is the
intrinsic object (`numerics/certify_lclm.py`, two independent proofs). **`Sym²(M)` has a rational
solution `R(t)`**, so:

> **What that clause does, explicitly.** "Removing the flat-band pole and using the reflection symmetry"
> is a short phrase for a construction, so here it is. With `x = 1/z` and `S(x) = Σₙ mₙ xⁿ = z·G(z)` the
> generating function of the moments themselves,
>
> ```
> S(x) = (1/3)/(1 + 2x) + Φ(x²/(1−x)²)/(1 − x),        Φ(t) = Σₘ νₘ tᵐ,
> ```
>
> exactly, checked coefficient by coefficient in exact rational arithmetic for `x⁰ … x²³⁰`. The first term
> is the flat-band pole; the gauge `1/(1−x) = ζ/x` is the Jacobian of writing `ζ·G_disp`; and
> **`t = x²/(1−x)² = 1/(z−1)²` is the degree-two invariant of the involution `x ↦ x/(2x−1)`, which is
> `E ↦ 2−E`**. So the passage to `t` is a **quadratic pullback, not a substitution**: every singular point
> of `M` has two preimages in the energy variable (`9t−1 → (4x−1)(2x+1)`, i.e. `E = 4` and `E = −2`;
> `5t−1 → 4x²+2x−1`; `4t−1 → (3x−1)(x+1)`; `t−1 → 2x−1`), and an annihilator guessed directly from the raw
> `mₙ` in the variable `x` has order **five**, not four. See
> [`numerics/verify_moment_bridge.py`](numerics/verify_moment_bridge.py) → `CERTIFICATE_moment_bridge.txt`.

- **Orthogonal Galois group.** The solution space of `M` carries a monodromy-invariant, nondegenerate,
  symmetric bilinear form; the differential Galois group is `G = O(3,ℂ)`, identity component
  `G° = SO(3,ℂ) ≅ PSL(2,ℂ)`. Equivalently `M` is homomorphic to its adjoint via an explicit order-2
  intertwiner. Hence **`M` *is* projectively equivalent to the symmetric square of a second-order
  operator `V₂`** (via a differential intertwiner — not a function multiplier or algebraic pullback, which
  is why the exponent test does not see it). The determinant character of the monodromy is the quadratic
  character of the genus-one twist curve `v² = (1−4t)(1−5t)(1−9t)`.
- **Modular parametrization — the closed form (proven).** `V₂` uniformizes `X(Γ₀(30)⁺)`:
  **`t = u/(u²+7u+1)`** with `u` the level-30 eta quotient above **generates the genus-zero function field
  of `X(Γ₀(30)⁺)`** (`1/t = u+7+1/u`; we call this a *Weber-like function parametrization*, reserving
  "Hauptmodul" for the `₂F₁` pullback `1728/j`), and the Schwarzian identity `{τ,t} = 2·Q_V(t)` holds exactly, with
  `Q_V = N(t)/[4t²(t−1)²(4t−1)²(5t−1)²(9t−1)²]`,
  `N = 24300t⁸−58860t⁷+73437t⁶−44294t⁵+15111t⁴−3160t³+407t²−30t+1`. The projective monodromy of `M` is
  the arithmetic lattice `Γ₀(30)⁺` (covolume `3π`, signature `(0;2,2,2,2,2;1 cusp)`); the five order-2
  points map to `t ∈ {1/9,1/5,1/4,1,∞}` and the cusp to `t=0`. ⇒ **hyperkagome LGF modular at level
  `30 = 2·3·5`**. Proven exactly in `numerics/certify_modular.py` → `CERTIFICATE_modular.txt`.
  The same uniformizing equation is tabulated as the row `Γ₀(30)⁺` of Lian–Yau (arXiv:hep-th/9507151)
  and as Conway–Norton class `30B` of Lian–Wiczer; `numerics/certify_tabulated.py` verifies both
  identifications exactly, Jacobian included. What is new here is the lattice-side identification and
  the pole-degree bound that makes the match a proof.
- **The weight-2 period `y₀ = Φ'/2` in explicit closed form.** `y₀ = [ρ₀(t)·W + ρ₁(t)·W'] / √((1−4t)(1−5t)(1−9t))`
  with `W = q·dt/dq` and explicit rational `ρ₀, ρ₁` — a **weight-two, depth-one quasimodular** form on `Γ₀(30)⁺`
  twisted by the determinant character. The derivative term is **provably essential**: `y₀` is *not* a
  (meromorphic modular form)×(algebraic function) at any weight or level, in particular **not an eta quotient**,
  so the period is genuinely quasimodular rather than modular. Certified exactly in `numerics/certify_y0.py` and
  `numerics/certify_y0_lemma.py` (→ `CERTIFICATE_y0.txt`, `CERTIFICATE_y0_lemma.txt`).
- **The generating function `Φ` lies outside the quasimodular module of `y₀` (exact obstruction).**
  Integrating the closed form once: `Φ = 2ρ₁·W/v + 2∫₀ᵗ Δ·(W/v) ds + 2/15` with
  **`Δ = ρ₀ − (ρ₁' − ρ₁·v'/v) = (13t−1)/(30t(t−1)²)`**. **Theorem:**
  `Φ ∉ ℚ(t) + v⁻¹(ℚ(t)W + ℚ(t)W' + ℚ(t)W'')` — no rational solution of `Ñ_v(S₂) = 2Δ` exists
  (indicial analysis forbids every pole and caps the degree at 3; the resulting 4-dim system is
  inconsistent, exactly over ℚ). The invariant content is the non-vanishing of the class
  `[2Δ] ∈ ℚ(t)/Ñ_v(ℚ(t))`; the individual representative `Δ` — and hence where its poles sit —
  depends on the choice of `v`-model and slice, so no structural meaning attaches to the particular
  form `(13t−1)/(30t(t−1)²)`. Physically: since `ζ·G_disp = Φ(t)` (`ζ = z−1`, `t = ζ⁻²`), the
  derivative `d[ζG_disp]/dz = −4y₀(t)/ζ³` **is** quasimodular, while `ζG_disp` itself is **not in
  that module**; the gap is an Eichler-type integral of the weight-4 character form `Δ·W²/v`.
  **Scope:** the theorem excludes the `ℚ(t)`-module displayed above. It does not by itself exclude
  quasimodularity with coefficients algebraic over `ℚ(t)`, nor the opposite `v`-parity.
  Certified in `numerics/certify_phi_obstruction.py` → `CERTIFICATE_phi_obstruction.txt`.
- **The ₂F₁ pullback is explicit and solvable in radicals (Atkin–Lehner Galois group).** The pullback
  `H = 1728/j` of the `₂F₁([1/12,5/12];[1];·)` representation has primitive minimal polynomial `P_H(H,t)` of
  bidegree **(8,72)** over `ℚ(t)`; its splitting field is the multiquadratic extension
  `ℚ(t)(√((1−t)(1−9t)), √((1−t)(1−5t)), √(1−4t))` — the full function field of `X₀(30)` — so
  `Gal(P_H) ≅ (ℤ/2)³` **is the Atkin–Lehner group of level 30**, and the pullback is solvable in radicals
  (explicit radical expression: M. van Hoeij, private communication, July 2026). The product of the three
  radicands is `((1−t)·v)²`: the determinant-character twist is the product of the quadratic layers.
  Proof-grade: annihilation is verified through `q^1200`, past the a-priori divisor bound `8·72+72·8 = 1152`,
  and the radical expression is checked as an **exact root** (polynomial identity, no series). For the
  `₂F₁([1/8,3/8];[1];·)` base the minimal polynomial has bidegree (4,24) on an index-two subfield.
  Data + certificate: `numerics/pullback_data.json`, `numerics/certify_pullback.py` → `CERTIFICATE_pullback.txt`.
- **Still non-Liouvillian ⇒ no algebraic or elementary closed form** — `G° = SO(3,ℂ)` is simple, hence
  non-solvable (an irreducible operator has Liouvillian solutions iff `G°` is solvable). This excludes
  algebraic/elementary forms only, and is fully consistent with the modular closed form (eta quotients and
  ₂F₁ are themselves non-Liouvillian).
- **The Watson reduction (concrete realization).** Exactly verified: `Φ(t) = 2/3 − (t/3) d/dt CT_θ log D₊` with
  `D₊ = (1−6t+3t²) − 2t²·Σcos2θᵢ − 8t^{3/2}·cosθ₁cosθ₂cosθ₃` — i.e. the hyperkagome LGF **is the classical
  generalized Watson integral** of the bcc lattice with 1st- (weight `t^{3/2}`) and 2nd- (sc shell, weight `t²`)
  neighbour hopping at spectral parameter `1−6t+3t²`. The `√t` double cover is exactly the `v²=(1−4t)(1−5t)(1−9t)`
  twist above; the singular set `{1/9,1/5,1/4,1}` are the band-map critical values. This anchors the level-30
  modular parametrization above on the lattice side (both describe the same object), and connects it to the
  classical second-neighbour cubic-lattice-LGF literature (Morita–Horiguchi, Glasser, Joyce 1998). Derived
  and verified in this repository (`numerics/verify_watson_reduction.py` → `CERTIFICATE_watson.txt`).

`M` is an order-3 Fuchsian period with six singular points (`t = 0, 1/9, 1/5, 1/4, 1, ∞`) plus an
apparent-type degree-7 locus. At `t = 0` (⇔ `z = ∞`) the local monodromy is **maximally unipotent (MUM)**:
exponents `{−1,−1,0}` (all integers), a single 3×3 Jordan block ⇒ maximal log power `n = 2`, and the Frobenius
basis takes the canonical MUM normal form `{y0, y0·log t + f1, y0·log²t/2 + f1·log t + f2}` after a rational
recombination (`numerics/verify_mum_normalform.py`). Note the exponents are
`{−1,−1,0}`, not the `{0,0,0}` of the Calabi–Yau normalization. It is **not** a ₃F₂ (which would have only three
singular points), and it is **not** of Calabi–Yau type (Zudilin: an order-3 CY operator would be a literal
symmetric square, and this is not). Headline exact special value: **`Re G(1) = 1/9`**.

## What is proven vs. what is open

- **Certified by exact computation** (rational arithmetic, no floating point in the decisive steps):
  `M` has no order-1 or order-2 right factor over ℚ(t); irreducibility over ℚ̄(t) then follows by a
  Galois-descent argument using the genuine log at `t = 0`; `M` is not a *literal* symmetric square
  (exponent triples are not arithmetic progressions); the Riemann scheme, including the exact degree-7
  exponents `{0,1,3}`.
- **Orthogonal structure, certified exactly** (`numerics/certify_orthogonal.py`,
  `CERTIFICATE_orthogonal.txt`): `Sym²(M)` has the rational solution `R(t)` (194 exact relations at two
  base points, made a proof by a Fuchs-relation budget with cap 109; Gram matrix nondegenerate); `M` is homomorphic to its adjoint via an explicit order-2
  intertwiner `T` (`rightremainder(M·T, adjoint(M)) = 0` exactly over ℚ(t)); the `t = 0` Jordan block is
  single/3×3 (`n = 2`); the determinant character is `−1` at `{1/9,1/5,1/4,∞}` ⇒ `G = O(3,ℂ)`, twist
  curve `v² = (1−4t)(1−5t)(1−9t)`. ⇒ `G° = SO(3,ℂ)`, projectively a symmetric square.
- **Where the `M ↔ Sym²(V₂)` intertwiner hides** (`numerics/certify_intertwiner.py`,
  `CERTIFICATE_intertwiner.txt`): the Wronskian log-derivative of `M` has residue `−5/2` at each of
  `t = 1/4, 1/5, 1/9`, so `Λ³(M)` is non-trivial while `Sym²(V₂)` (projective normal form) is unimodular ⇒
  **no rational homomorphism exists in either direction**, necessarily. The single conjugation
  `M_v = v·M·v⁻¹` shifts that log-derivative by `3v′/v` and kills the order-two character, after which the
  intertwiner is the **order-one** `T = ρ₀ + ρ₁ d/dt` (remainder of `M_v·T` mod `Sym²(V₂)` exactly `0`).
  The same parity (`det Sym² = (det)⁴`) removes the obstruction between `Sym²(M)` and `Sym⁴(V₂)`, whose
  homomorphy over ℚ(t) is easy to check in Maple and is **not** certified here.
- **Modular parametrization, certified exactly** (`numerics/certify_modular.py`,
  `CERTIFICATE_modular.txt`): the level-30 eta quotient `u` satisfies the Ligozat conditions (⇒ modular
  function on `Γ₀(30)`); `t = u/(u²+7u+1)` matches the `V₂` MUM mirror map to 80 orders; and the Schwarzian
  identity `{τ,t} = 2·Q_V(t)` holds — an a priori degree-≤8 rational identity verified far past its bound,
  hence a proof. ⇒ `V₂` uniformizes `X(Γ₀(30)⁺)`; `t` generates its function field; LGF modular at level 30.
- **Verified to the guess-and-verify standard of the field:** `M` annihilates the symmetry-reduced
  series exactly over ℚ through `t^111` — **112 relations against 57 free coefficients, a margin of
  55, zero residuals** — plus an independent from-scratch reconstruction and an independent
  Bloch-Hamiltonian moment computation.
- **Open (acknowledged):** an unconditional creative-telescoping *proof* that `M` is the minimal
  annihilator of the LGF (see [`numerics/CT_SETUP.md`](numerics/CT_SETUP.md)). The weight-2 period `y₀`
  is now given in explicit quasimodular closed form (above), and is **proven not** to be an eta quotient or a
  (modular form)×(algebraic function) at any weight — so that is no longer an open item.

## Reproduce

Requires only Python 3.10+ with NumPy, SciPy, SymPy, mpmath (no SageMath):

```bash
pip install -r requirements.txt
python numerics/verify.py                    # moment fingerprint + exact annihilation of M (margin 55)
python numerics/certify_factor.py            # no order-1/2 factor over Q; not-literal-Sym^2; irreducibility
python numerics/certify_lclm.py              # L4 = M d/dt is the UNIQUE factorization (no LCLM splitting)
python numerics/certify_orthogonal.py        # Sym^2(M) rational solution => G = O(3,C); intertwiner; n=2; det char
python numerics/certify_intertwiner.py       # det character is the ONLY obstruction to a rational M <-> Sym^2(V2)
python numerics/certify_modular.py           # V2 uniformizes X(Gamma_0(30)+): t=u/(u^2+7u+1) generates the field; Schwarzian
python numerics/certify_tabulated.py         # that uniformizing equation IS the tabulated Gamma_0(30)+ / class-30B entry
python numerics/certify_nonliouvillian.py    # genuine log at t=0 => non-Liouvillian (no algebraic/elementary form)
python numerics/certify_p7_apparent.py       # p7 is an APPARENT locus (all 3 local solutions log-free)
python numerics/certify_bridge.py            # explicit V2 + conic point; bridge f0^2 = P(y0), V2(f0)=0 to t^107
python numerics/certify_y0.py                # y0 = Phi'/2 closed form (weight-2 depth-1 quasimodular)
python numerics/certify_y0_lemma.py          # y0 is NOT (modular form of any weight) x (algebraic)
python numerics/certify_phi_obstruction.py   # Phi lies outside y0's quasimodular module; Delta identity + theorem
python numerics/certify_pullback.py          # 2F1 pullback 1728/j: minimal polynomial (8,72); solvable in radicals;
                                             #   Galois group = Atkin-Lehner (Z/2)^3; proof past the divisor bound 1152
python numerics/verify_mum_normalform.py     # t=0 is MUM: canonical normal form, n=2
python numerics/verify_watson_reduction.py   # LGF = generalized bcc(1,2) Watson integral
python numerics/strengthen_certification.py  # overdetermination margin + independent Bloch moments
python numerics/verify_specialvalues.py      # Re G(1) = 1/9 by symmetric BZ quadrature
python numerics/verify_vanhove_log.py        # log-divergent van Hove point t=1 (E=0,2) vs smooth controls
python numerics/vm_crosscheck.py             # matches the Varma–Monien spectrum & 1/(t_VM+1) pole
python numerics/verify_moment_bridge.py      # raw moments m_n <-> Phi: S(x) display, quadratic pullback t=x^2/(1-x)^2
```

All twenty-one exit `0` on system Python.

Each script validates its primitives on operators of known structure, with negative controls,
before the certified computation. Five of them (`certify_bridge`, `certify_intertwiner`,
`certify_lclm`, `certify_phi_obstruction`, `certify_pullback`) write their certificate file
themselves; every other `numerics/CERTIFICATE*.txt` is the captured standard output of the
script of the same name, e.g.

```bash
python numerics/certify_modular.py > numerics/CERTIFICATE_modular.txt
```

## Layout

```
paper/           main.tex, main.pdf, figs/         — the manuscript (journal class)
paper_arxiv/     main.tex, main.pdf, make_arxiv.py — the same paper on the plain article class;
                                                     make_arxiv.py --check asserts the two bodies are
                                                     byte-identical
numerics/
  lattice.pkl                                       — hyperkagome unit cell (12 sites, directed NN bonds)
  moments230.json, nu.json                          — exact integer moments m0..m230; symmetry-reduced nu
  M_coeffs.json                                      — the certified order-3 operator M (integer coeffs)
  V2_data.json                                       — the explicit order-2 operator V2 (p, q, Q_V), the
                                                       transported N = D^3+B2 D^2+B1 D+B0, and the conic point
  pullback_data.json                                 — the (8,72) minimal polynomial of the 2F1 pullback 1728/j,
                                                       van Hoeij's radical expression, and the (4,24) equation
                                                       for the [1/8,3/8] base
  verify.py, certify_factor.py, certify_lclm.py, certify_orthogonal.py, certify_intertwiner.py,
  certify_modular.py,
  certify_nonliouvillian.py, certify_p7_apparent.py, certify_bridge.py, certify_y0.py, certify_y0_lemma.py,
  certify_phi_obstruction.py, certify_pullback.py, certify_tabulated.py, verify_mum_normalform.py,
  verify_watson_reduction.py,
  strengthen_certification.py, verify_specialvalues.py, verify_vanhove_log.py, vm_crosscheck.py,
  verify_moment_bridge.py
                                                     — reproduction / certification scripts (21, all exit 0)
                                                       (certify_modular.py: level-30 uniformization;
                                                        certify_tabulated.py: the same equation as the entry
                                                        tabulated by Lian–Yau and by Lian–Wiczer;
                                                        certify_bridge.py: explicit V2 + the bridge identity)
  extend_moments.py                                  — closed-walk moment generator (provenance)
  CERTIFICATE*.txt                                   — generated certificates
  CT_SETUP.md, CERTIFY.md                            — creative-telescoping route; DFactor/Magma cross-check
  M_maple.txt, M_magma.txt, M_sage.py                — M in other CAS syntaxes (for independent factorization)
```

## Methodology

The operator was obtained by differential approximants (the Guttmann-style "guess from exact moments,
then verify" method), and all certifications were carried out by computer algebra in exact rational
arithmetic with SymPy. Because the certification is by exact verification, it does not depend on how
`M` was originally guessed. Every reported result is reproduced by the scripts above.

## License

Code and data: MIT ([`LICENSE`](LICENSE)). The manuscript text and figures in `paper/` are © the
authors, all rights reserved.

## Citation

See [`CITATION.cff`](CITATION.cff).
