# Hyperkagome lattice Green's function: modular uniformization at level 30

This repository accompanies the manuscript *"Lattice Green's function of the hyperkagome lattice:
modular uniformization at level 30 from an orthogonal differential Galois group"*
([`paper/main.pdf`](paper/main.pdf)). It contains the exact operator, the data it was built from,
and self-contained scripts that reproduce every certification.

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
> (i.e. `1/t = u+7+1/u`) is a **Hauptmodul of `Γ₀(30)⁺`**, and `{τ,t} = 2·Q_V(t)` holds exactly. So the
> hyperkagome LGF is **modular at level `30 = 2·3·5`** — apparently the first lattice Green's function at a
> modular level with three distinct prime factors. Proven exactly in
> [`numerics/certify_modular.py`](numerics/certify_modular.py) → `CERTIFICATE_modular.txt`.

## Result

Varma & Monien (*Lattice Green's functions for kagome, diced and hyperkagome lattices*,
[Phys. Rev. E **87**, 032109 (2013)](https://doi.org/10.1103/PhysRevE.87.032109); arXiv:1211.5666)
reduced the 3D hyperkagome (Na₄Ir₃O₈) density of states to a threefold integral (their Eq. 16) and
wrote that they "found no way to exactly solve" it.

Working from exact lattice moments, we determine that the hyperkagome Green's function — after
removing the flat-band pole and using the exact reflection symmetry about `E = 1` — is annihilated
by an **irreducible order-3 Picard–Fuchs operator `M` of degree 15**, and that `M` is **not a
*literal* symmetric square** of a second-order operator (its Frobenius exponent triples are not the
arithmetic progressions a symmetric square would force). **Nevertheless `Sym²(M)` has a rational
solution `R(t)`**, so:

- **Orthogonal Galois group.** The solution space of `M` carries a monodromy-invariant, nondegenerate,
  symmetric bilinear form; the differential Galois group is `G = O(3,ℂ)`, identity component
  `G° = SO(3,ℂ) ≅ PSL(2,ℂ)`. Equivalently `M` is homomorphic to its adjoint via an explicit order-2
  intertwiner. Hence **`M` *is* projectively equivalent to the symmetric square of a second-order
  operator `V₂`** (via a differential intertwiner — not a function multiplier or algebraic pullback, which
  is why the exponent test does not see it). The determinant character of the monodromy is the quadratic
  character of the genus-one twist curve `v² = (1−4t)(1−5t)(1−9t)`.
- **Modular parametrization — the closed form (proven).** `V₂` uniformizes `X(Γ₀(30)⁺)`:
  **`t = u/(u²+7u+1)`** with `u` the level-30 eta quotient above is a **Hauptmodul of `Γ₀(30)⁺`**
  (`1/t = u+7+1/u`), and the Schwarzian identity `{τ,t} = 2·Q_V(t)` holds exactly, with
  `Q_V = N(t)/[4t²(t−1)²(4t−1)²(5t−1)²(9t−1)²]`,
  `N = 24300t⁸−58860t⁷+73437t⁶−44294t⁵+15111t⁴−3160t³+407t²−30t+1`. The projective monodromy of `M` is
  the arithmetic lattice `Γ₀(30)⁺` (covolume `3π`, signature `(0;2,2,2,2,2;1 cusp)`); the five order-2
  points map to `t ∈ {1/9,1/5,1/4,1,∞}` and the cusp to `t=0`. ⇒ **hyperkagome LGF modular at level
  `30 = 2·3·5`**. Proven exactly in `numerics/certify_modular.py` → `CERTIFICATE_modular.txt`.
- **The weight-2 period `y₀ = Φ'/2` in explicit closed form.** `y₀ = [ρ₀(t)·W + ρ₁(t)·W'] / √((1−4t)(1−5t)(1−9t))`
  with `W = q·dt/dq` and explicit rational `ρ₀, ρ₁` — a **weight-two, depth-one quasimodular** form on `Γ₀(30)⁺`
  twisted by the determinant character. The derivative term is **provably essential**: `y₀` is *not* a
  (meromorphic modular form)×(algebraic function) at any weight or level, in particular **not an eta quotient**,
  so the period is genuinely quasimodular rather than modular. Certified exactly in `numerics/certify_y0.py` and
  `numerics/certify_y0_lemma.py` (→ `CERTIFICATE_y0.txt`, `CERTIFICATE_y0_lemma.txt`).
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
basis takes the canonical MUM normal form `{y0, y0·log t + f1, y0·log²t/2 + f1·log t + f2}` after an integer
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
- **Modular parametrization, certified exactly** (`numerics/certify_modular.py`,
  `CERTIFICATE_modular.txt`): the level-30 eta quotient `u` satisfies the Ligozat conditions (⇒ modular
  function on `Γ₀(30)`); `t = u/(u²+7u+1)` matches the `V₂` MUM mirror map to 80 orders; and the Schwarzian
  identity `{τ,t} = 2·Q_V(t)` holds — an a priori degree-≤8 rational identity verified far past its bound,
  hence a proof. ⇒ `V₂` uniformizes `X(Γ₀(30)⁺)`; `t` is a Hauptmodul; LGF modular at level 30.
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
python numerics/certify_orthogonal.py        # Sym^2(M) rational solution => G = O(3,C); intertwiner; n=2; det char
python numerics/certify_modular.py           # V2 uniformizes X(Gamma_0(30)+): t=u/(u^2+7u+1) Hauptmodul; Schwarzian
python numerics/certify_nonliouvillian.py    # genuine log at t=0 => non-Liouvillian (no algebraic/elementary form)
python numerics/certify_p7_apparent.py       # p7 is an APPARENT locus (all 3 local solutions log-free)
python numerics/strengthen_certification.py  # overdetermination margin + independent Bloch moments
python numerics/verify_specialvalues.py      # Re G(1) = 1/9 by symmetric BZ quadrature
python numerics/vm_crosscheck.py             # matches the Varma–Monien spectrum & 1/(t_VM+1) pole
```

Each `certify_*` script validates its primitives on operators with known structure before the real
run, and writes a plain-text certificate (`numerics/CERTIFICATE*.txt`).

## Layout

```
paper/           main.tex, main.pdf, figs/         — the manuscript
numerics/
  lattice.pkl                                       — hyperkagome unit cell (12 sites, directed NN bonds)
  moments230.json, nu.json                          — exact integer moments m0..m230; symmetry-reduced nu
  M_coeffs.json                                      — the certified order-3 operator M (integer coeffs)
  verify.py, certify_factor.py, certify_orthogonal.py, certify_modular.py, certify_nonliouvillian.py,
  certify_p7_apparent.py, strengthen_certification.py, verify_specialvalues.py, vm_crosscheck.py
                                                     — reproduction / certification scripts
                                                       (certify_modular.py: level-30 Hauptmodul / uniformization)
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
authors, all rights reserved (pending journal submission).

## Citation

See [`CITATION.cff`](CITATION.cff).
