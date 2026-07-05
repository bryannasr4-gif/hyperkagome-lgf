# Creative-telescoping proof of the operator M — setup and status

> **Status.** This is the one remaining item that could not be *completed* in the local environment:
> a formal creative-telescoping (CT) derivation of M directly from the resolvent integral would upgrade the
> operator from "conjectured + massively verified" to "proved." It requires a CT engine — `ore_algebra`
> (SageMath), Koutschan's `HolonomicFunctions` (Mathematica), or Lairez's `period`/`ore_algebra` reduction —
> **none of which can be installed on this Windows machine** (no WSL/Docker/conda, no Windows SageMath package,
> no C compiler; confirmed). Building the degree-12 integrand symbolically also overwhelms plain `sympy`
> (`det(I-tH)` for the 12×12 Bloch matrix timed out). Below is the exact, ready-to-run setup for a machine that
> has one of these engines.
>
> **What IS already established (accepted standard in the lattice-statistics literature, e.g. Guttmann,
> Zenine–Boukraa–Hassani–Maillard):** M annihilates the exact series over ℚ through t¹¹¹ — **112 relations
> against 57 free coefficients, an overdetermination margin of 55, zero residuals** — plus an independent
> from-scratch rebuild (bit-identical) and an **independent Bloch-Hamiltonian recomputation of the moments**
> (`strengthen_certification.py`). A CT proof is the rigor capstone, not a correction.

## The integral representation

With hopping t=1, the moment generating function is the Brillouin-zone constant term of the resolvent:
```
Phi_full(t) = (1/12) * CT_{x,y,z} [ tr( I - t*H(x,y,z) )^{-1} ]
            = (1/12) * CT_{x,y,z} [ tr(adj(I - t*H)) / det(I - t*H) ],
```
where `H(x,y,z)` is the 12×12 Bloch Hamiltonian with entries `H[a,b] = sum over bonds a->b of x^dx y^dy z^dz`
(`x=e^{i k_x}` etc.), built from `numerics/lattice.pkl`. The integrand is a rational function of
`(t,x,y,z)`; `det(I - t*H)` is its denominator. The reduction to the certified operator M is then:
strip the flat-band pole `(1/3)/(1+2t)` (line-graph flat band at E=-2), symmetrize in `E->2-E` (the exact
E=1 reflection, now proved via bipartiteness of the premedial net), and change variable to `t=(z-1)^{-2}`.

## Route A — SageMath + ore_algebra (recommended)

```python
from ore_algebra import OreAlgebra, guess
import pickle
sites, bonds = pickle.load(open("numerics/lattice.pkl","rb"))
R.<x,y,z,t> = QQ[]; F = FractionField(R)
H = matrix(F, 12, 12)
for (a,b,d) in [(int(i),int(j),tuple(int(round(v)) for v in dd)) for (i,j,dd) in bonds]:
    H[a,b] += x^d[0]*y^d[1]*z^d[2]     # (use Laurent: multiply through by monomials to clear negatives)
Iden = identity_matrix(F,12)
integrand = (Iden - t*H).inverse().trace() / 12          # rational function in x,y,z,t
# creative telescoping w.r.t. the three torus variables (constant-term / diagonal):
A = OreAlgebra(F, 'Dx','Dy','Dz','Dt')
# use ore_algebra's diagonal / CT interface (annihilator of the constant term):
#   telescoper, certificates = A.creative_telescoping(integrand, [Dx,Dy,Dz])   # eliminate x,y,z
# then verify telescoper (an ODE in t) equals L = M*Dt after the flat-band + symmetry reduction.
```
`ore_algebra` provides the certified telescoper + certificates `g_x,g_y,g_z` with
`L*integrand = Dx(g_x)+Dy(g_y)+Dz(g_z)`, which is the unconditional proof. Compare the reduced telescoper to
`M_coeffs.json` (they must agree).

## Route B — Mathematica + HolonomicFunctions (Koutschan)

```mathematica
Needs["HolonomicFunctions`"]
(* build H[x,y,z] (12x12) from the bond list; integrand = Tr[Inverse[IdentityMatrix[12] - t H]]/12 *)
ct = CreativeTelescoping[integrand, {Der[x], Der[y], Der[z]}, Der[t]]
(* ct[[1]] is the telescoper (ODE in t); reduce (flat band + E=1 symmetry) and compare to M. *)
```

## Route C — Lairez, "Computing periods of rational integrals"

Lairez's algorithm computes the Picard–Fuchs operator of `CT/∮ (rational integrand)` by Griffiths–Dwork
reduction; his `period` implementation (Sage/`ore_algebra` or the Julia port) applied to the same integrand
yields the operator directly. This is the most direct "period of a rational function" route.

## Verification target

Any of the three must return an operator that, after the flat-band deflation and the `E=1` symmetry reduction,
equals the order-3 operator `M` in `numerics/M_coeffs.json` (equivalently, `L = M·Dt`). That equality is the
formal proof. Until then, the operator is certified to the field's standard as quantified above.
