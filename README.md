# Hyperkagome lattice Green's function: an irreducible order-3 Picard–Fuchs operator

This repository accompanies the manuscript *"Lattice Green's function of the hyperkagome lattice:
an irreducible third-order Picard–Fuchs operator and the absence of an elliptic closed form"*
([`paper/main.pdf`](paper/main.pdf)). It contains the exact operator, the data it was built from,
and self-contained scripts that reproduce every certification.

## Result

Varma & Monien (*Lattice Green's functions for kagome, diced and hyperkagome lattices*,
[Phys. Rev. E **87**, 032109 (2013)](https://doi.org/10.1103/PhysRevE.87.032109); arXiv:1211.5666)
reduced the 3D hyperkagome (Na₄Ir₃O₈) density of states to a threefold integral (their Eq. 16) and
wrote that they "found no way to exactly solve" it.

Working from exact lattice moments, we determine that the hyperkagome Green's function — after
removing the flat-band pole and using the exact reflection symmetry about `E = 1` — is annihilated
by an **irreducible order-3 Picard–Fuchs operator `M` of degree 15**, and that `M` is **not a
symmetric square** of a second-order (elliptic) operator. Two independent consequences follow:

- **Not a symmetric square ⇒ no closed form in complete elliptic integrals** (unlike the SC/BCC/FCC
  cubic-lattice LGFs, which *are* symmetric squares of elliptic operators).
- **A genuine logarithmic solution at `t = 0` ⇒ `M` is non-Liouvillian ⇒ no algebraic or elementary
  closed form.** (Non-Liouvillian alone does *not* exclude an elliptic form — `K` itself is
  non-Liouvillian — so the elliptic exclusion comes solely from the not-symmetric-square property.
  The two legs together give: no closed form in algebraic, elementary, or elliptic functions.)

`M` is an order-3 Fuchsian period with six singular points (`t = 0, 1/9, 1/5, 1/4, 1, ∞`) plus an
apparent-type degree-7 locus, of Calabi–Yau type in the broad sense — **not** a ₃F₂ (which would have
only three singular points). Headline exact special value: **`Re G(1) = 1/9`**.

## What is proven vs. what is open

- **Certified by exact computation** (rational arithmetic, no floating point in the decisive steps):
  `M` has no order-1 or order-2 right factor over ℚ(t); irreducibility over ℚ̄(t) then follows by a
  Galois-descent argument using the genuine log at `t = 0`; `M` is not a symmetric square at every
  singular point; the Riemann scheme, including the exact degree-7 exponents `{0,1,3}`.
- **Verified to the guess-and-verify standard of the field:** `M` annihilates the symmetry-reduced
  series exactly over ℚ through `t^111` — **112 relations against 57 free coefficients, a margin of
  55, zero residuals** — plus an independent from-scratch reconstruction and an independent
  Bloch-Hamiltonian moment computation.
- **Open (acknowledged):** an unconditional creative-telescoping *proof* that `M` is the minimal
  annihilator of the LGF (the whole no-closed-form conclusion is conditional on this). See
  [`numerics/CT_SETUP.md`](numerics/CT_SETUP.md).

## Reproduce

Requires only Python 3.10+ with NumPy, SciPy, SymPy, mpmath (no SageMath):

```bash
pip install -r requirements.txt
python numerics/verify.py                    # moment fingerprint + exact annihilation of M (margin 55)
python numerics/certify_factor.py            # no order-1/2 factor over Q; not-Sym^2; irreducibility
python numerics/certify_nonliouvillian.py    # genuine log at t=0 => non-Liouvillian
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
  verify.py, certify_factor.py, certify_nonliouvillian.py, certify_p7_apparent.py,
  strengthen_certification.py, verify_specialvalues.py, vm_crosscheck.py
                                                     — reproduction / certification scripts
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
author, all rights reserved (pending journal submission).

## Citation

See [`CITATION.cff`](CITATION.cff).
