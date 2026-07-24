"""certify_bridge.py -- the explicit second-order operator V2 and the bridge identity.

This certificate makes public the object on which Secs. V-VII of the paper rest: the
explicit order-two operator V2 over Q(t) whose symmetric square is projectively M, the
rational conic point that produced it, and the *bridge identity* that ties V2 back to the
lattice Green's function in a single exact computation.

WHAT IS PROVED HERE (all in exact rational arithmetic; no floats, no SageMath):

  (a) N := D^3 + B2 D^2 + B1 D + B0, the transported symmetric-square operator, is a
      LITERAL symmetric square: writing p = B2/3 and q = (B1 - 2p^2 - p')/4, the
      symmetric-square consistency condition  B0 == 4pq + 2q'  holds identically.
      Hence V2 = D^2 + p D + q is defined over Q(t).

  (b) Q_V of Eq. (9) of the paper is the projective normal form of that V2:
      Q_V == q - p'/2 - p^2/4, identically.  (Q_V is gauge invariant, so this is the
      statement that the paper's Eq. (9) and the V2 emitted here are the same operator
      up to a function multiplier.)

  (c) THE BRIDGE.  With P = v0 + v1 D + v2 D^2 the order-two intertwiner built from the
      rational conic point, and y0 = Phi'/2 the lattice period (read from the certified
      moments in nu.json, NOT from any V2 data), the function

              f0 := sqrt( P(y0) / P(y0)|_{t=0} )

      satisfies V2(f0) = 0 exactly, as a power series, to the full length supported by
      the moment file.  This single identity validates the conic point, the intertwiner
      and V2 simultaneously, and it is what connects the modular Secs. VI-VII to the
      lattice.

VALIDATION FIRST (per the project's verify-before-run rule): every series primitive used
below is exercised on operators of known structure before the real run --
  * the symmetric-square criterion is checked to ACCEPT a manufactured Sym^2 and to
    REJECT a nearby non-symmetric-square,
  * the series square root is checked against a known square,
  * the ODE residual routine is checked to return 0 on an exact solution of a control
    operator and nonzero on a perturbation of it.

Exit code 0 and 'ALL CHECKS PASS' means every item above verified.
"""
import json
import os
from fractions import Fraction as F

import sympy as sp

t = sp.symbols('t')
HERE = os.path.dirname(os.path.abspath(__file__))
LINES = []


def say(s=""):
    print(s)
    LINES.append(s)


def check(label, ok):
    say("  %-68s %s" % (label, "PASS" if ok else "FAIL"))
    if not ok:
        raise SystemExit("FAILED: " + label)
    return ok


# ----------------------------------------------------------------- primitives
def poly_coeffs(e):
    P = sp.Poly(sp.expand(e), t)
    if P.is_zero:
        return [F(0)]
    out = [F(0)] * (P.degree() + 1)
    for (k,), c in P.as_dict().items():
        out[k] = F(int(sp.numer(c)), int(sp.denom(c)))
    return out


def deriv(a):
    return [F(n + 1) * a[n + 1] for n in range(len(a) - 1)]


def sermul(a, b, N):
    r = [F(0)] * N
    for i, ai in enumerate(a):
        if ai and i < N:
            for j, bj in enumerate(b):
                if bj and i + j < N:
                    r[i + j] += ai * bj
    return r


def sersqrt(g, N):
    """principal square root of a unit series (g[0] == 1)."""
    assert g[0] == 1
    h = [F(1)] + [F(0)] * (N - 1)
    for n in range(1, N):
        h[n] = (g[n] - sum(h[k] * h[n - k] for k in range(1, n))) / 2
    return h


def laurent(expr, N):
    """Laurent expansion of a rational function at t=0: (coeffs, valuation)."""
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    nc, dc = poly_coeffs(num), poly_coeffs(den)
    vd = next(i for i, c in enumerate(dc) if c != 0)
    dc = dc[vd:]
    inv = [F(0)] * N
    inv[0] = 1 / dc[0]
    for n in range(1, N):
        inv[n] = -sum(dc[k] * inv[n - k]
                      for k in range(1, min(n, len(dc) - 1) + 1)) / dc[0]
    return sermul(nc, inv, N), -vd


def ode_residual(pp, qq, h, N):
    """Laurent coefficients of h'' + p h' + q h, as a dict {power: coeff}."""
    hp, hpp = deriv(h), deriv(deriv(h))
    pser, pval = laurent(pp, N + 4)
    qser, qval = laurent(qq, N + 4)
    E = {}

    def dep(ser, val):
        for i, c in enumerate(ser):
            if c:
                E[i + val] = E.get(i + val, F(0)) + c

    dep(hpp, 0)
    for i, c in enumerate(sermul(pser, hp, N)):
        if c:
            E[i + pval] = E.get(i + pval, F(0)) + c
    for i, c in enumerate(sermul(qser, h, N)):
        if c:
            E[i + qval] = E.get(i + qval, F(0)) + c
    return E


def symsq_defect(B2e, B1e, B0e):
    """0 iff D^3+B2 D^2+B1 D+B0 is a literal symmetric square."""
    pp = sp.cancel(B2e / 3)
    qq = sp.cancel((B1e - 2 * pp ** 2 - sp.diff(pp, t)) / 4)
    return sp.cancel(sp.together(4 * pp * qq + 2 * sp.diff(qq, t) - B0e)), pp, qq


# ------------------------------------------------------- 0. validate primitives
say("=" * 78)
say("certify_bridge.py -- explicit V2 and the bridge identity to the lattice period")
say("=" * 78)
say()
say("[0] VALIDATION OF PRIMITIVES ON KNOWN STRUCTURE (before any real computation)")

# (0a) symmetric-square criterion: build Sym^2 of a known V2 = D^2 + pD + q.
p_c = sp.cancel(1 / (t - 3))
q_c = sp.cancel((5 * t + 2) / (t ** 2 + 1))
B2_c = 3 * p_c
B1_c = 2 * p_c ** 2 + sp.diff(p_c, t) + 4 * q_c
B0_c = 4 * p_c * q_c + 2 * sp.diff(q_c, t)
d_c, p_r, q_r = symsq_defect(B2_c, B1_c, B0_c)
check("(0a) criterion ACCEPTS a manufactured symmetric square", d_c == 0)
check("(0a) and recovers p, q exactly",
      sp.cancel(p_r - p_c) == 0 and sp.cancel(q_r - q_c) == 0)
d_bad, _, _ = symsq_defect(B2_c, B1_c, B0_c + t)
check("(0a) criterion REJECTS a perturbed (non-Sym^2) operator", d_bad != 0)

# (0b) series square root
base = [F(1), F(3), F(-2), F(7), F(1), F(0), F(0), F(0)]
sq = sermul(base, base, 8)
root = sersqrt(sq, 8)
check("(0b) series sqrt recovers a known square", root == base)

# (0c) ODE residual on the Euler control y'' - (1/t) y' + (1/t^2) y = 0, exact solution y = t
ctrl = [F(0), F(1)] + [F(0)] * 20
res = ode_residual(sp.cancel(-1 / t), sp.cancel(1 / t ** 2), ctrl, 18)
check("(0c) residual == 0 on an exact solution of the control operator",
      all(c == 0 for c in res.values()))
res_bad = ode_residual(sp.cancel(-1 / t), sp.cancel(1 / t ** 2 + 1), ctrl, 18)
check("(0c) residual != 0 on a perturbed control operator",
      any(c != 0 for c in res_bad.values()))
say()

# ------------------------------------------------------------- 1. load V2 data
D = json.load(open(os.path.join(HERE, "V2_data.json")))


def rat(name):
    return sp.cancel(sp.sympify(D[name]["num"]) / sp.sympify(D[name]["den"]))


p, q, Qpub = rat("p"), rat("q"), rat("Q")
B2, B1, B0 = rat("B2"), rat("B1"), rat("B0")
v0, v1, v2 = rat("v0"), rat("v1"), rat("v2")

say("[1] THE DATA  (numerics/V2_data.json)")
say("    N = D^3 + B2 D^2 + B1 D + B0: the order-three operator produced by the conic")
say("    transport of the invariant quadratic form, INDEPENDENTLY of p and q below.")
say("    deg numer/denom:  B2 = %d/%d,  B1 = %d/%d,  B0 = %d/%d"
    % (sp.Poly(sp.numer(B2), t).degree(), sp.Poly(sp.denom(B2), t).degree(),
       sp.Poly(sp.numer(B1), t).degree(), sp.Poly(sp.denom(B1), t).degree(),
       sp.Poly(sp.numer(B0), t).degree(), sp.Poly(sp.denom(B0), t).degree()))
say("    candidate V2 = D^2 + p D + q:  p = %d/%d,   q = %d/%d"
    % (sp.Poly(sp.numer(p), t).degree(), sp.Poly(sp.denom(p), t).degree(),
       sp.Poly(sp.numer(q), t).degree(), sp.Poly(sp.denom(q), t).degree()))
say("    conic point P = v0 + v1 D + v2 D^2, deg v0,v1,v2 = %d,%d,%d"
    % (sp.Poly(v0, t).degree(), sp.Poly(v1, t).degree(), sp.Poly(v2, t).degree()))
say()

# --------------------------------------- 2. (a) literal symmetric square, (b) Q_V
say("[2] SYMMETRIC-SQUARE CONSISTENCY AND THE PROJECTIVE NORMAL FORM")
say("    (B2, B1, B0 are taken from the transported operator in the data file; the")
say("    criterion below is a genuine constraint on that trio, not a tautology.)")
defect, p_back, q_back = symsq_defect(B2, B1, B0)
check("(a) N satisfies the literal-Sym^2 criterion B0 == 4pq + 2q'", defect == 0)
check("(a) the p, q recovered from N match the stored V2 exactly",
      sp.cancel(p_back - p) == 0 and sp.cancel(q_back - q) == 0)

Qmine = sp.cancel(q - sp.diff(p, t) / 2 - p ** 2 / 4)
check("(b) Q_V == q - p'/2 - p^2/4  (projective normal form)",
      sp.cancel(sp.together(Qmine - Qpub)) == 0)

# independent tie to the paper's Eq. (9), embedded here as literal integers
NUMER = [1, -30, 407, -3160, 15111, -44294, 73437, -58860, 24300]
paperQ = sp.cancel(sum(c * t ** i for i, c in enumerate(NUMER))
                   / (4 * t ** 2 * (t - 1) ** 2 * (4 * t - 1) ** 2
                      * (5 * t - 1) ** 2 * (9 * t - 1) ** 2))
check("(b) Q_V equals Eq. (9) of the paper as printed",
      sp.cancel(sp.together(Qpub - paperQ)) == 0)
say()

# --------------------------------------------------------------- 3. the bridge
say("[3] THE BRIDGE:  f0 = sqrt(P(y0)/P(y0)|_0)  solves V2(f0) = 0")
nu = [F(a, b) for a, b in json.load(open(os.path.join(HERE, "nu.json")))]
y0 = [F(n + 1) * nu[n + 1] / 2 for n in range(len(nu) - 1)]   # y0 = Phi'/2
say("    y0 = Phi'/2 read from nu.json: %d exact terms, head %s"
    % (len(y0), [str(x) for x in y0[:4]]))

y0p, y0pp = deriv(y0), deriv(deriv(y0))
NS = len(y0) - 2
g = [F(0)] * NS
for pol, ser in ((v0, y0), (v1, y0p), (v2, y0pp)):
    gm = sermul(poly_coeffs(pol), ser, NS)
    g = [g[i] + gm[i] for i in range(NS)]
check("    P(y0) is a unit at t=0 (so the square root is a power series)", g[0] != 0)

gn = [x / g[0] for x in g]
f0 = sersqrt(gn, NS)
check("    the extracted square root squares back exactly",
      sermul(f0, f0, NS) == gn)

E = ode_residual(p, q, f0, NS - 2)
worst = max((k for k, c in E.items() if c != 0), default=None)
nz = sum(1 for c in E.values() if c != 0)
ORD_OK = NS - 6
tested = [k for k in E if k <= ORD_OK]
check("    V2(f0) == 0 on every Laurent coefficient through order t^%d "
      "(%d coefficients)" % (ORD_OK, len(tested)),
      all(E[k] == 0 for k in tested))
say("    (nonzero residual coefficients anywhere in the window: %d%s)"
    % (nz, "" if worst is None else "; highest index %d = truncation edge" % worst))
say()

# ------------------------------------------------------------- 4. negative control
say("[4] NEGATIVE CONTROL: the same test must FAIL for a perturbed conic point")
gbad = [F(0)] * NS
for pol, ser in ((v0 + t ** 3, y0), (v1, y0p), (v2, y0pp)):
    gm = sermul(poly_coeffs(pol), ser, NS)
    gbad = [gbad[i] + gm[i] for i in range(NS)]
fbad = sersqrt([x / gbad[0] for x in gbad], NS)
Eb = ode_residual(p, q, fbad, NS - 2)
check("    perturbing v0 by t^3 breaks the bridge",
      any(Eb.get(k, F(0)) != 0 for k in range(0, ORD_OK)))
say()

say("=" * 78)
say("ALL CHECKS PASS")
say("=" * 78)

with open(os.path.join(HERE, "CERTIFICATE_bridge.txt"), "w") as f:
    f.write("\n".join(LINES) + "\n")
print("\nwrote CERTIFICATE_bridge.txt")
