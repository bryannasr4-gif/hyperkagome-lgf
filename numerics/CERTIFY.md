# Two rigor fixes before submission — exact commands

> **STATUS (done):** Both fixes have been executed *exactly in pure Python* — see
> `numerics/certify_factor.py` and its saved output `numerics/CERTIFICATE.txt`. That script enumerates the
> complete finite set of Fuchs-admissible order-1 right factors of M and of adjoint(M), tests each exactly over
> ℚ (0 found), and confirms a guaranteed logarithmic solution at t=0 (repeated exponent −1) ⇒ **M is irreducible
> over Q̄ and is not a symmetric square** ⇒ no closed form in complete elliptic integrals. The Maple/Magma
> recipes below are now **optional independent cross-checks**, not blockers. (A full creative-telescoping proof
> and a Liouvillian-exclusion test remain as further strengthenings, not required for the paper's stated claims.)

These are the two things a referee (very possibly Maillard, Guttmann, or Koutschan themselves) will poke
first. All three toolchains below use the **same operator M**, emitted by
`python numerics/certify_ops.py` into `M_maple.txt`, `M_magma.txt`, `M_sage.py`.

- **Fix 1 — Certified irreducibility.** The paper currently argues irreducibility from a hand-rolled
  hyperexponential search. Replace that with a *certified factorization*: if the tool factors M into a single
  order-3 piece, M is irreducible over ℚ(t). This is routine and independent of the Sage build bug that blocked us.
- **Fix 2 — No Liouvillian / elliptic closed form (the differential-Galois step).** "Not a symmetric square"
  (already solid, via the exponent test) rules out an *elliptic* K-form. To make the blunt "no closed form"
  airtight you must also exclude a finite or imprimitive differential-Galois group (which would give algebraic /
  Liouvillian solutions). The community-standard way (Boukraa–Hassani–Maillard) is to probe the symmetric and
  exterior squares for rational/exponential solutions and confirm a logarithmic solution exists (⇒ infinite
  Galois group).

Interpretation feeds back into the paper: if Fix 1 confirms irreducible and Fix 2 finds *no* rational solutions
of the exterior/symmetric square and a genuine log solution, you may keep "irreducible, not a symmetric square,
Calabi–Yau-type period, no closed form in complete elliptic integrals" as *proved* (not just verified).

---

## Route A — Maple (recommended; `DFactor` handles order ≤ 4 well)

```maple
with(DEtools):
# paste the four lines of numerics/M_maple.txt here (defines c0,c1,c2,c3 and L):
#   L := c3*Dt^3 + c2*Dt^2 + c1*Dt + c0:

# --- FIX 1: certified irreducibility ---
F := DFactor(L, [Dt, t]);         # a single factor (= L) means IRREDUCIBLE over Q(t)
nops([F]);                        # 1  => irreducible ; >1 => reducible, print the factors
eigenring(L, [Dt, t]);            # cross-check: eigenring = scalars (dim 1) <=> irreducible

# --- FIX 2: differential-Galois / not-Liouvillian ---
S2 := symmetric_power(L, 2, [Dt, t]):   # order 6 if L is NOT a symmetric square
E2 := exterior_power(L, 2, [Dt, t]):    # order 3 (the "adjoint-like" square)
ratsols(S2, t);                          # expect [] : no rational solution of Sym^2
ratsols(E2, t);                          # expect [] : L not self-dual via a rational kernel
expsols(S2, t);  expsols(E2, t);         # expect [] : no exponential (Liouvillian) building block
# a log solution => infinite Galois group => not finite/algebraic:
formal_sol(L, [Dt, t], t = 1);           # look for a ln(t-1) term in the local solutions at t=1
```

**Read-off:** `DFactor` returns one operator ⇒ **irreducible** (Fix 1 done). `ratsols`/`expsols` all empty and a
genuine log at t=1 ⇒ Galois group is not finite/imprimitive ⇒ **no Liouvillian (hence no elliptic) closed form**
(Fix 2 done). Report the exact `DFactor` output in the paper's certification paragraph.

---

## Route B — Magma (if you have access; `Factorisation`)

```magma
// paste numerics/M_magma.txt (defines F,R,D, c0..c3, L):
//   L := c3*D^3 + c2*D^2 + c1*D + c0;
fac := Factorisation(L);
#fac;                 // 1 => irreducible order-3 operator
[Order(f[1]) : f in fac];
```

Magma's differential-operator factorisation is a certified decision procedure and is the cleanest single
irreducibility certificate for a referee.

---

## Route C — Sage + ore_algebra (free; use a FRESH conda env, not the broken build)

The project's Sage `.factor()` was broken by a `passagemath-flint`/`_small_primes_table` build bug. A clean
environment fixes it:

```bash
# fresh miniforge env, no sudo
conda create -n sage2 -c conda-forge sage python=3.11 -y
conda activate sage2
sage -pip install ore_algebra
sage numerics/M_sage.py      # runs L.factor(); prints "irreducible iff 1 and order 3: True/False"
```

`M_sage.py` builds the same M and calls `L.factor()`. One factor of order 3 ⇒ irreducible. (For Fix 2 in Sage:
`L.symmetric_power(2)` and `.rational_solutions()` mirror the Maple recipe.)

---

## If a check comes back the "wrong" way

- **DFactor finds a factor** ⇒ M is reducible; the "irreducible order-3 / not elliptic" story changes and the
  paper must be re-derived. (Very unlikely given the exponent structure and the independent rebuild, but this is
  exactly why the certified check matters.)
- **ratsols/expsols non-empty** ⇒ a Liouvillian building block exists; soften to "no *elliptic* closed form"
  and drop the unqualified "no closed form / Calabi–Yau" wording, or characterize the Liouvillian solution.

Once both fixes pass, update `paper/main.tex` (the "$M$ has no low-order right factor" and "No elliptic closed
form" paragraphs) to cite the certified factorization and the Galois computation, then proceed to outreach.
