"""
CERTIFY: the differential Galois group of the hyperkagome operator M is ORTHOGONAL.

This certificate CORRECTS and REPLACES the paper's earlier negative claim ("M is not a
symmetric square in any equivalent sense => no closed form in complete elliptic integrals;
SL3-type Galois group"). The correction is due to J.-M. Maillard (private communication,
July 2026): although M is irreducible and is NOT a *literal* symmetric square, the symmetric
square Sym^2(M) has a RATIONAL solution, so the solution space of M carries a monodromy-
invariant nondegenerate quadratic form and the differential Galois group G is contained in the
orthogonal group O(3,C).  Equivalently M is homomorphic to its adjoint by an order-2 intertwiner.
Hence M IS projectively equivalent to the symmetric square of a second-order operator, and an
elliptic/2F1 closed form is EXPECTED rather than excluded (the explicit second-order operator
lives on the genus-one curve u^2 = (1-4t)(1-5t)(1-9t)).

  ---------------------------------------------------------------------------------------------
  ATTRIBUTION.  The three facts this script certifies -- (i) Homomorphisms(adjoint(M),M) != 0
  with the explicit order-2 intertwiner T, (ii) the rational solution R of Sym^2(M), and
  (iii) the resulting O(3,C) structure -- are due to J.-M. Maillard (private communication,
  July 2026).  This script is an INDEPENDENT exact re-verification of his results.
  ---------------------------------------------------------------------------------------------

All decisive steps are exact (rational arithmetic / operator algebra over Q(t); no floating
point).  Every routine is validated on operators with KNOWN structure before the real run.

  (A) Sym^2 MEMBERSHIP, exact, at TWO ordinary base points t0 = 1/2 and t0 = -1/3:
      R = sum c_ij y_i y_j with CONSTANT c_ij (the y_i a fundamental system of M).  Solved
      from 6 coefficients, verified on 194 further exact relations at each base point.
      This is a PROOF and not merely a margin, by a FUCHS-RELATION BUDGET: if
      F = R - sum c_ij y_i y_j were not identically zero, the span W of {y_i y_j, R} would be
      7-dimensional and its minimal operator Fuchsian of order 7; every local exponent of that
      operator is bounded below by the certified data (pairwise sums of M's exponents at each
      singular point, the pole/zero orders of R, and DISTINCT nonnegative integers at every
      apparent point), and the Fuchs relation then caps the vanishing order of any nonzero
      element of W at an M-ordinary point at V <= 109.  The verified vanishing to order 200
      exceeds the budget, so F == 0 identically.  (The budget is recomputed and printed by
      this script.)  The Gram matrix C = (c_ij) is symmetric and nondegenerate (rank 3).

  (B) INTERTWINER, fully symbolic (exact over Q(t)):  with Maillard's order-2 operator T,
      rightremainder(M o T, adjoint(M)) = 0 exactly.  So T maps solutions of adjoint(M) to
      solutions of M; a nonzero homomorphism between the irreducible M and adjoint(M) is an
      isomorphism, i.e. M is self-dual.  Self-dual + irreducible + ODD order => the invariant
      form is symmetric and nondegenerate => G is contained in O(3,C).  (This upgrades the
      margin of (A) to an operator-level proof and is independent of it.)

  (C) JORDAN STRUCTURE at t=0 (= z=infinity; a MUM point -- see verify_mum_normalform.py):
      indicial polynomial
      -64*rho*(rho+1)^2 => exponents {-1,-1,0}, all INTEGERS => local monodromy is unipotent;
      the number of log-free local solutions is exactly 1 = the number of Jordan blocks, so the
      unipotent has a single 3x3 block and the maximal power of log is n = 2 (log^2 present).
      This settles the formerly-open n=1-vs-2 item: n = 2 (EVEN), consistent with the orthogonal
      case of the Hassani-Maillard-Zenine parity conjecture (arXiv:2502.05543).

  (D) DETERMINANT CHARACTER (which orthogonal group):  the Wronskian log-derivative -c2/c3 has
      residue (sum of local exponents) - 3 at each singular point; this is a HALF-integer
      exactly at t = 1/9, 1/5, 1/4, infinity, so det(monodromy) = -1 there and +1 at 0, 1, p7 --
      an EVEN number (four) of det=-1 points, matching the branch set of
      u^2 = (1-4t)(1-5t)(1-9t).  So G = O(3,C) in full (not just SO(3,C)); the determinant
      character is the quadratic character of that genus-one curve.  Combined with (C) (a genuine
      log => G infinite and non-semisimple), the identity component is G^0 = SO(3,C).

  (E) INPUTS to the Clifford / Lie-Kolchin argument that pins G^0 = SO(3,C) (written out in the
      paper) are the ALREADY-CERTIFIED facts: irreducibility over Qbar(t) (certify_factor.py /
      CERTIFICATE.txt) and the genuine t=0 logarithm (certify_nonliouvillian.py, re-counted here
      in (C)).  Non-Liouvillianity now has the SHARPER proof: G^0 = SO(3,C) is simple, hence
      non-solvable, hence M has no Liouvillian solutions (still no algebraic/elementary form).

Run:  python certify_orthogonal.py          (writes the certificate to stdout)
"""
import json
import os
import sys
from collections import defaultdict
from fractions import Fraction as F

import sympy as sp

t = sp.symbols('t')
HERE = os.path.dirname(os.path.abspath(__file__))

# =====================================================================================
# Exact ascending-coefficient polynomial / power-series toolkit (Fraction arithmetic)
# =====================================================================================
def pmul(a, b):
    r = [F(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    r[i + j] += ai * bj
    return r

def pscale(a, s):
    return [ai * s for ai in a]

def pshift(p, a):
    """coefficients of p(s + a), exact Taylor shift, truncated to deg p."""
    n = len(p)
    out = [F(0)] * n
    for c in reversed(p):
        new = [F(0)] * (len(out) + 1)
        for i, oi in enumerate(out):
            new[i + 1] += oi
            new[i] += oi * a
        new[0] += c
        out = new[:n]
    return out

def pdesc(coeffs_desc):
    return [F(c) for c in reversed(coeffs_desc)]

def series_inv(d, N):
    assert d[0] != 0
    inv = [F(0)] * N
    inv[0] = 1 / d[0]
    for n in range(1, N):
        s = F(0)
        for k in range(1, min(n, len(d) - 1) + 1):
            s += d[k] * inv[n - k]
        inv[n] = -s / d[0]
    return inv

def series_mul(a, b, N):
    r = [F(0)] * N
    for i in range(min(len(a), N)):
        if a[i]:
            for j in range(min(len(b), N - i)):
                if b[j]:
                    r[i + j] += a[i] * b[j]
    return r

def ff(n, i):
    r = F(1)
    for k in range(i):
        r *= (n - k)
    return r

def gauss_solve(A, b):
    """Exact solve of a square Fraction system; returns None if singular."""
    n = len(A)
    Mx = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = next((r for r in range(col, n) if Mx[r][col] != 0), None)
        if piv is None:
            return None
        Mx[col], Mx[piv] = Mx[piv], Mx[col]
        inv = Mx[col][col]
        Mx[col] = [v / inv for v in Mx[col]]
        for r in range(n):
            if r != col and Mx[r][col] != 0:
                f = Mx[r][col]
                Mx[r] = [Mx[r][k] - f * Mx[col][k] for k in range(n + 1)]
    return [Mx[i][n] for i in range(n)]

# =====================================================================================
# Ore-algebra primitives (operator = ascending list [c0,c1,...] meaning sum ci D^i)
# =====================================================================================
def trim(op):
    op = [sp.cancel(c) for c in op]
    while len(op) > 1 and op[-1] == 0:
        op.pop()
    return op

def Dcomp(B):
    out = [sp.Integer(0)] * (len(B) + 1)
    for k, b in enumerate(B):
        out[k] += sp.diff(b, t)
        out[k + 1] += b
    return trim(out)

def opmul(A, B):
    res = [sp.Integer(0)]
    DiB = [c for c in B]
    for i, a in enumerate(A):
        if a != 0:
            term = [sp.cancel(a * c) for c in DiB]
            if len(term) > len(res):
                res = res + [sp.Integer(0)] * (len(term) - len(res))
            for k, c in enumerate(term):
                res[k] += c
        DiB = Dcomp(DiB)
    return trim(res)

def opright_divmod(A, B):
    A = [sp.cancel(c) for c in A]
    m = len(B) - 1
    bm = B[m]
    Q = [sp.Integer(0)]
    while len(A) - 1 >= m and any(c != 0 for c in A):
        n = len(A) - 1
        an = A[n]
        q = sp.cancel(an / bm)
        deg = n - m
        if len(Q) <= deg:
            Q = Q + [sp.Integer(0)] * (deg + 1 - len(Q))
        Q[deg] += q
        qD = [sp.Integer(0)] * deg + [q]
        sub = opmul(qD, B)
        newlen = max(len(A), len(sub))
        A = A + [sp.Integer(0)] * (newlen - len(A))
        sub = sub + [sp.Integer(0)] * (newlen - len(sub))
        A = [sp.cancel(A[k] - sub[k]) for k in range(newlen)]
        while len(A) > 1 and A[-1] == 0:
            A.pop()
    return trim(Q), trim(A)

def right_remainder(A, B):
    return opright_divmod(A, B)[1]

def adjoint(L):
    n = len(L) - 1
    out = [sp.Integer(0)] * (n + 1)
    for i in range(n + 1):
        ai = L[i]
        for k in range(i + 1):
            out[k] += (-1)**i * sp.binomial(i, k) * sp.diff(ai, t, i - k)
    return trim([sp.cancel(c) for c in out])

def sym2_order2(L2):
    """Classical symmetric square of an order-2 operator (returns order-3 operator)."""
    a0, a1 = sp.cancel(L2[0] / L2[2]), sp.cancel(L2[1] / L2[2])
    return [sp.cancel(4 * a0 * a1 + 2 * sp.diff(a0, t)),
            sp.cancel(2 * a1**2 + sp.diff(a1, t) + 4 * a0),
            sp.cancel(3 * a1),
            sp.Integer(1)]

def indicial_exponents(op, p):
    """Local exponents of `op` at the finite regular singular point t=p (exact)."""
    ssym = sp.symbols('sig')
    r = sp.symbols('r')
    terms = defaultdict(lambda: sp.Integer(0))
    for i, c in enumerate(op):
        ci = sp.cancel(c).subs(t, p + ssym)
        if ci == 0:
            continue
        coeff, m = ci.as_leading_term(ssym).as_coeff_exponent(ssym)
        fall = sp.prod([r - a for a in range(i)]) if i > 0 else sp.Integer(1)
        terms[sp.nsimplify(m) - i] += coeff * fall
    mn = min(terms.keys())
    roots = sp.roots(sp.Poly(sp.expand(terms[mn]), r))
    exps = []
    for rt, mult in roots.items():
        exps += [sp.nsimplify(rt)] * mult
    return sorted(exps, key=lambda x: sp.re(sp.N(x)))

def exponents_at_infinity(op):
    """Exponents at t=infinity in the paper convention (exponent rho <=> y ~ t^{-rho})."""
    s = sp.symbols('s')
    n = len(op) - 1
    expr = sum(sp.expand(c) * sp.ff(-s, i) * t**(n - i) for i, c in enumerate(op))
    poly = sp.Poly(sp.expand(expr), t)
    ind = poly.coeff_monomial(t**poly.degree())
    roots = sp.roots(sp.Poly(sp.expand(ind), s))
    exps = []
    for rt, mult in roots.items():
        exps += [sp.nsimplify(-rt)] * mult
    return sorted(exps, key=lambda z: sp.re(sp.N(z)))

def count_log_free(op, p, K=8):
    """(number of log-free formal local solutions at t=p, sorted exponents)."""
    x = sp.symbols('x')
    s = sp.symbols('s')
    n = len(op) - 1
    cij = []
    for i in range(n + 1):
        d = {}
        for (j,), co in sp.Poly(sp.expand(sp.expand(op[i]).subs(t, p + x)), x).terms():
            if co != 0:
                d[j] = sp.Rational(co)
        cij.append(d)
    qd = defaultdict(lambda: sp.Integer(0))
    for i in range(n + 1):
        for j, co in cij[i].items():
            qd[j - i] += co * sp.prod([s - a for a in range(i)])
    dmin = min(dd for dd in qd if sp.expand(qd[dd]) != 0)
    roots = sp.roots(sp.Poly(sp.expand(qd[dmin]), s))
    exps = []
    for rt, m in roots.items():
        exps += [sp.nsimplify(rt)] * m
    emin = min(exps, key=lambda z: sp.re(sp.N(z)))
    rows = []
    for r in range(K):
        row = []
        for k in range(K):
            dd = dmin + r - k
            val = qd[dd].subs(s, emin + k) if (k <= r and dd in qd) else sp.Integer(0)
            row.append(sp.Rational(sp.expand(val)))
        rows.append(row)
    nullity = K - sp.Matrix(rows).rank()
    return nullity, sorted(exps, key=lambda z: sp.re(sp.N(z)))

# =====================================================================================
# Load M  (both as Fraction coefficient lists and as a sympy operator)
# =====================================================================================
Mj = json.load(open(os.path.join(HERE, "M_coeffs.json")))
assert Mj["order"] == 3 and len(Mj["coeffs"]) == 4
cF = [[F(int(x)) for x in lst] for lst in Mj["coeffs"]]                 # ascending Fractions
Msp = [sum(int(c) * t**j for j, c in enumerate(Mj["coeffs"][i])) for i in range(4)]

xm1 = [F(-1), F(1)]

# =====================================================================================
# Validate the primitives on synthetic operators with KNOWN structure
# =====================================================================================
print("validating operator/series primitives on known operators...")
# adjoint involutive
_L = [t**2 + 1, t, sp.Integer(3), sp.Integer(1)]
assert all(sp.cancel(a - b) == 0 for a, b in zip(_L, adjoint(adjoint(_L)))), "adjoint not involutive"
# classical sym2 sanity
assert trim(sym2_order2([sp.Integer(-1), sp.Integer(0), sp.Integer(1)])) == \
       trim([sp.Integer(0), -sp.Integer(4), sp.Integer(0), sp.Integer(1)]), "sym2(D^2-1) != D^3-4D"
# right division: exact factor => zero remainder
_B = [t + 1, sp.Integer(1)]
_T = opmul([sp.Integer(1), t, sp.Integer(1)], _B)
assert right_remainder(_T, _B) == [sp.Integer(0)], "right-division on a known factor failed"
# log-free counter: t^3 D^3 has 3 (apparent), theta^3 has 1 (single Jordan block, n=2)
assert count_log_free([sp.Integer(0), sp.Integer(0), sp.Integer(0), t**3], sp.Integer(0))[0] == 3
assert count_log_free([sp.Integer(0), t, 3 * t**2, t**3], sp.Integer(0))[0] == 1
print("  [validate] adjoint / sym2 / right-division / log-counter: PASS\n")

R_str = "-(1/272)*(15t^2+17t-8)^2 / ( t^2 (t-1)^2 (4t-1)(5t-1)(9t-1) )"

# =====================================================================================
# PART A:  Sym^2 membership -- R = sum c_ij y_i y_j, exact, margin >= 180, TWO base points
# =====================================================================================
print("=" * 78)
print("(A) Sym^2 MEMBERSHIP  (R is a rational solution of Sym^2(M))")
print("    R(t) =", R_str)
print("=" * 78)

def sym2_membership(t0, N):
    cs = [pshift(ci, t0) for ci in cF]
    assert cs[3][0] != 0, "base point must be ordinary (c3(t0) != 0)"

    def solve_series(ic):
        a = [F(0)] * N
        a[0], a[1], a[2] = ic
        for m in range(N - 3):
            s = F(0)
            for i in range(4):
                for j, bij in enumerate(cs[i]):
                    if bij:
                        idx = m - j + i
                        if 0 <= idx <= m + 2:
                            s += bij * ff(idx, i) * a[idx]
            a[m + 3] = -s / (cs[3][0] * ff(m + 3, 3))
        return a

    Y = [solve_series(ic) for ic in ([F(1), F(0), F(0)],
                                     [F(0), F(1), F(0)],
                                     [F(0), F(0), F(1)])]
    num = pscale(pmul(pdesc([15, 17, -8]), pdesc([15, 17, -8])), F(-1, 272))
    den = pmul(pmul(pmul(pmul([F(0), F(0), F(1)], pmul(xm1, xm1)),
                         [F(-1), F(4)]), [F(-1), F(5)]), [F(-1), F(9)])
    R = series_mul(pshift(num, t0), series_inv(pshift(den, t0), N), N)
    pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
    P = [series_mul(Y[i], Y[j], N) for (i, j) in pairs]
    A6 = [[P[k][m] for k in range(6)] for m in range(6)]
    b6 = [R[m] for m in range(6)]
    csol = gauss_solve(A6, b6)
    if csol is None:
        return None
    bad = [m for m in range(N) if sum(csol[k] * P[k][m] for k in range(6)) != R[m]]
    C = [[csol[0], csol[1] / 2, csol[2] / 2],
         [csol[1] / 2, csol[3], csol[4] / 2],
         [csol[2] / 2, csol[4] / 2, csol[5]]]
    det = (C[0][0] * (C[1][1] * C[2][2] - C[1][2] * C[2][1])
           - C[0][1] * (C[1][0] * C[2][2] - C[1][2] * C[2][0])
           + C[0][2] * (C[1][0] * C[2][1] - C[1][1] * C[2][0]))
    return csol, bad, det

N = 200
brief_c = [F(-17, 84), F(-460, 441), F(44927, 3087),
           F(-1925200, 157437), F(83066560, 1102059), F(-712021969, 7714413)]
okA = True
for t0 in [F(1, 2), F(-1, 3)]:
    res = sym2_membership(t0, N)
    assert res is not None, "leading 6x6 system singular at t0=%s" % t0
    csol, bad, det = res
    margin = N - 6
    good = (len(bad) == 0 and det != 0 and margin >= 180)
    okA = okA and good
    print("  base point t0 = %-5s : margin = %d exact relations, residuals bad = %s" %
          (t0, margin, bad[:3] if bad else "none"))
    print("       c=(c11,c12,c13,c22,c23,c33) =", [str(v) for v in csol])
    print("       Gram det(C) = %s   rank = %d   =>  %s" %
          (det, 3 if det != 0 else "<3",
           "nondegenerate symmetric form (monodromy in O(3,C))" if det != 0 else "DEGENERATE"))
    if t0 == F(1, 2):
        match = (csol == brief_c and det == F(3437476900, 7714413))
        print("       matches Maillard/BRIEF reference values at t0=1/2:", "YES" if match else "NO")
        okA = okA and match
# ---- Fuchs-relation budget: why the margin-194 agreement is a PROOF ------------------
# If F = R - sum c_ij y_i y_j were nonzero, W = span{y_i y_j, R} would be 7-dimensional and
# its minimal operator L7 Fuchsian of order 7.  Fuchs relation (n=7, n(n-1)/2 = 21):
#     sum_{finite sing s} (expsum_{L7}(s) - 21)  +  (expsum_{L7}(inf) + 21)  =  0.
# The exponent multiset of the minimal operator of a function space dominates (sorted,
# entrywise) the valuation multiset of any spanning family, so each term is bounded below by
#   * pairwise sums of M's certified exponents + val(R), at each singular point of M;
#   * 7 DISTINCT nonnegative integers at any apparent point (all solutions holomorphic there):
#     excess >= 0, and at the base point t0, where F vanishes to order V, sum >= 15 + V;
#   * pairwise sums of {3/2,2,3} plus 3 (R ~ t^-3), at infinity (z = 1/t convention).
exps_fin = {sp.Integer(0): indicial_exponents(Msp, sp.Integer(0)),
            sp.Rational(1, 9): indicial_exponents(Msp, sp.Rational(1, 9)),
            sp.Rational(1, 5): indicial_exponents(Msp, sp.Rational(1, 5)),
            sp.Rational(1, 4): indicial_exponents(Msp, sp.Rational(1, 4)),
            sp.Integer(1): indicial_exponents(Msp, sp.Integer(1))}
valR = {sp.Integer(0): -2, sp.Rational(1, 9): -1, sp.Rational(1, 5): -1,
        sp.Rational(1, 4): -1, sp.Integer(1): -2}          # pole orders of the factored R
budget = sp.Integer(0)
for p, ex in exps_fin.items():
    pair = sum(ex[i] + ex[j] for i in range(3) for j in range(i, 3))
    budget += (pair + valR[p]) - 21
budget += 7 * (sp.Integer(0 + 1 + 2 + 3 + 4 + 6) - 21)     # p7: exponents {0,1,3} (certified
                                                           # log-free in certify_p7_apparent.py),
                                                           # pairwise sums {0,1,2,3,4,6}; val(R)=0
inf_ex = sorted([-e for e in exponents_at_infinity(Msp)],   # helper returns the negated values;
                key=lambda z: sp.re(sp.N(z)))               # paper convention: y ~ t^{-rho}
assert inf_ex == [sp.Rational(3, 2), sp.Integer(2), sp.Integer(3)], inf_ex
pair_inf = sum(inf_ex[i] + inf_ex[j] for i in range(3) for j in range(i, 3))   # = 26
budget += (pair_inf + 3) + 21                               # R ~ t^{-3} at infinity
Vmax = 6 - budget                                           # 0 >= (15+V-21) + budget => V <= 6-budget
print("  Fuchs-relation budget: if F = R - sum c_ij y_i y_j were nonzero, the order-7 minimal")
print("  operator of span{y_i y_j, R} would violate the Fuchs relation unless the vanishing")
print("  order V of F at the (M-ordinary) base point obeyed V <= %s." % Vmax)
print("  Verified vanishing order >= %d > %s  =>  F == 0 identically: the margin is a PROOF." % (N, Vmax))
okA = okA and (N > Vmax)
print("  => R is a rational solution of Sym^2(M); the invariant form C is symmetric,")
print("     nondegenerate (rank 3) at BOTH base points.   (A):", "PASS" if okA else "FAIL")

# =====================================================================================
# PART B:  intertwiner remainder identity  (fully symbolic, exact over Q(t))
# =====================================================================================
print()
print("=" * 78)
print("(B) INTERTWINER  Homomorphisms(adjoint(M), M) != 0  (M is self-dual)")
print("=" * 78)

Madj = adjoint(Msp)

# ---- validate the remainder machinery on a KNOWN self-dual operator (positive control) ----
Msc = [sp.Integer(0), -sp.Integer(4), sp.Integer(0), sp.Integer(1)]     # D^3 - 4D = sym2(D^2-1)
assert right_remainder(opmul(Msc, [sp.Integer(1)]), adjoint(Msc)) == [sp.Integer(0)], \
    "positive control: identity is a known intertwiner of the self-dual D^3-4D"
print("  [validate] D^3-4D self-dual, T0=1 is an intertwiner (rem=0): PASS")

# ---- Maillard's explicit order-2 intertwiner T (email, private communication July 2026) ----
x = t
p7 = 101025*x**7 - 369600*x**6 + 455798*x**5 - 290956*x**4 + 93657*x**3 - 17580*x**2 + 1800*x - 64
quartic = 225*x**4 + 510*x**3 + 49*x**2 - 272*x + 64
a2 = sp.Rational(1, 22730625) * (x - 1) * p7 * quartic
a1 = (22*x**11 - sp.Rational(50528, 1347)*x**10 - sp.Rational(2363078, 33675)*x**9
      + sp.Rational(243633628, 1515375)*x**8 - sp.Rational(809086106, 22730625)*x**7
      - sp.Rational(3189283384, 22730625)*x**6 + sp.Rational(3695323154, 22730625)*x**5
      - sp.Rational(1919500508, 22730625)*x**4 + sp.Rational(185285356, 7576875)*x**3
      - sp.Rational(93296624, 22730625)*x**2 + sp.Rational(8706304, 22730625)*x
      - sp.Rational(349184, 22730625))
a0 = sp.Rational(2, 7576875) * (416728125*x**11 - 841083750*x**10 - 1197586350*x**9
      + 3977462330*x**8 - 2088650740*x**7 - 2674062330*x**6 + 4392973398*x**5
      - 2776823266*x**4 + 957803671*x**3 - 187304696*x**2 + 19545656*x - 868288) / (x - 1)
T = [sp.cancel(a0), sp.cancel(a1), sp.cancel(a2)]

MT = opmul(Msp, T)
rem = right_remainder(MT, Madj)
okB = all(sp.cancel(sp.together(c)) == 0 for c in rem)
print("  order(M o T) = %d, order(adjoint(M)) = %d, remainder order <= %d" %
      (len(MT) - 1, len(Madj) - 1, len(rem) - 1))
print("  rightremainder(M o T, adjoint(M)) == 0 exactly over Q(t):", "PASS" if okB else "FAIL")
# non-vacuity: perturbing T by a nonzero constant breaks the identity
rem_pert = right_remainder(opmul(Msp, [T[0] + 1, T[1], T[2]]), Madj)
nonvac = any(sp.cancel(c) != 0 for c in rem_pert)
print("  [control] perturbed T (a0+1) gives NONZERO remainder:", "PASS" if nonvac else "FAIL")
print("  => T: solutions(adjoint(M)) -> solutions(M) is a nonzero homomorphism; M irreducible")
print("     => it is an isomorphism => M is self-dual => invariant bilinear form on Sol(M).")
print("     ODD order + irreducible => the form is SYMMETRIC & nondegenerate => G in O(3,C).")
okB = okB and nonvac

# =====================================================================================
# PART C:  Jordan structure at t=0 (a MUM point; see verify_mum_normalform.py) -- n = 2
# =====================================================================================
print()
print("=" * 78)
print("(C) JORDAN STRUCTURE at t=0  (t=0 <=> z=infinity; MUM -- single maximal unipotent block)")
print("=" * 78)
rho = sp.symbols('rho')
lo1, lo2, lo3 = cF[1][0], cF[2][1], cF[3][2]     # lowest coeffs contributing to the indicial poly
ind = sp.expand(sp.Integer(int(lo3)) * rho * (rho - 1) * (rho - 2)
                + sp.Integer(int(lo2)) * rho * (rho - 1) + sp.Integer(int(lo1)) * rho)
ind_roots = sp.roots(sp.Poly(ind, rho))
nlf, exps0 = count_log_free(Msp, sp.Integer(0))
all_integer = all(sp.nsimplify(e).is_integer for e in exps0)
n_log = 3 - nlf                                   # single block of size 3 => max log power 2
okC = (sp.factor(ind) == sp.factor(-64 * rho * (rho + 1)**2)
       and ind_roots == {sp.Integer(-1): 2, sp.Integer(0): 1}
       and nlf == 1 and all_integer and n_log == 2)
print("  indicial polynomial at t=0 :", sp.factor(ind), " (= -64*rho*(rho+1)^2)")
print("  exponents {-1,-1,0} :", [str(e) for e in exps0], " all integers =>",
      "unipotent local monodromy" if all_integer else "NOT unipotent")
print("  log-free local solutions = %d of 3  => number of Jordan blocks = %d" % (nlf, nlf))
print("  [control] theta^3 (single 3x3 block, n=2) log-free =",
      count_log_free([sp.Integer(0), t, 3*t**2, t**3], sp.Integer(0))[0], "of 3")
print("  => a single 3x3 Jordan block => maximal power of log is n = %d (log^2 present)." % n_log)
print("     Even n = 2 is consistent with the ORTHOGONAL case of arXiv:2502.05543.   (C):",
      "PASS" if okC else "FAIL")

# =====================================================================================
# PART D:  determinant character -- which orthogonal group (O(3) vs SO(3)) + twist curve
# =====================================================================================
print()
print("=" * 78)
print("(D) DETERMINANT CHARACTER  (Wronskian residues => det of monodromy => G = O(3,C))")
print("=" * 78)
wlog = sp.cancel(-Msp[2] / Msp[3])                # W'/W  for the Wronskian W of M
fin_sings = [sp.Integer(0), sp.Rational(1, 9), sp.Rational(1, 5), sp.Rational(1, 4), sp.Integer(1)]
det_minus = []
okD = True
print("  point      res(-c2/c3)   sum(exponents)-3   det(monodromy)")
for p in fin_sings:
    res = sp.residue(wlog, t, p)
    ex = indicial_exponents(Msp, p)
    w_exp = sum(ex) - 3
    match = (res == w_exp)
    okD = okD and match
    half = (sp.nsimplify(2 * res) % 2 != 0)
    det = -1 if half else 1
    if det == -1:
        det_minus.append(p)
    print("  t=%-7s  %-11s   %-16s   %+d%s" %
          (p, res, w_exp, det, "" if match else "   [MISMATCH]"))
# p7 locus: exponents {0,1,3} => sum-3 = 1 (integer) => det = +1 (also apparent)
print("  p7 locus   (exps {0,1,3}) sum-3 = 1 (integer)         +1   (apparent)")
# infinity: exponents {3/2,2,3} (paper convention) => sum = 13/2 half-integer => det -1
inf_exps = exponents_at_infinity(Msp)
sum_inf = sum(inf_exps)
inf_half = (sp.nsimplify(2 * sum_inf) % 2 != 0)
if inf_half:
    det_minus.append(sp.oo)
print("  t=inf      (exps %s) sum = %s  => det = %s" %
      ([str(e) for e in inf_exps], sum_inf, "-1" if inf_half else "+1"))
# cross-check: residues over ALL of P^1 sum to zero
u = sp.symbols('u')
res_inf = sp.residue((-Msp[2] / Msp[3]).subs(t, 1 / u) * (-1 / u**2), u, 0)
print("  [cross-check] residue at infinity of (-c2/c3)dt = %s (half-integer: %s)" %
      (res_inf, sp.nsimplify(2 * res_inf) % 2 != 0))
branch = {sp.Rational(1, 9), sp.Rational(1, 5), sp.Rational(1, 4), sp.oo}
det_minus_set = set(det_minus)
okD = okD and (det_minus_set == branch) and (len(det_minus_set) % 2 == 0)
print("  det(monodromy) = -1 exactly at", sorted([str(s) for s in det_minus_set]),
      "(%d points, EVEN)" % len(det_minus_set))
print("  => G is NOT inside SO(3,C); G = O(3,C) in full.  The determinant character is the")
print("     quadratic character of the genus-one curve  u^2 = (1-4t)(1-5t)(1-9t)")
print("     (branch set {1/4,1/5,1/9,infinity} = the det=-1 set).   (D):", "PASS" if okD else "FAIL")

# =====================================================================================
# PART E:  inputs to the Clifford / Lie-Kolchin argument (already-certified elsewhere)
# =====================================================================================
print()
print("=" * 78)
print("(E) INPUTS to the Clifford/Lie-Kolchin pin-down of G^0 = SO(3,C) (proof in the paper)")
print("=" * 78)
print("  irreducible over Qbar(t)                 : certify_factor.py / CERTIFICATE.txt")
print("  genuine logarithm at t=0 (=> G infinite, : certify_nonliouvillian.py; re-counted in")
print("     non-semisimple, so G^0 non-solvable)     (C) above (log-free dim 1 of 3)")
print("  invariant form symmetric & nondegenerate : (A),(B) above (Schur: unique up to scale;")
print("                                              odd dim => symmetric, not symplectic)")
print("  => the only non-solvable connected subgroup of O(3,C) is SO(3,C), so G^0 = SO(3,C)")
print("     (~= PSL(2,C) via Sym^2 of the standard rep).  Non-Liouvillian now follows from")
print("     G^0 = SO(3,C) simple => non-solvable (sharper than the old not-Sym^2 route).")

# =====================================================================================
# CERTIFICATE
# =====================================================================================
allpass = okA and okB and okC and okD
print()
print("================ ORTHOGONAL / SO(3) CERTIFICATE ================")
print("(A) Sym^2(M) has the rational solution R (margin >=%d, 2 base points, C nondeg): %s"
      % (N - 6, "PASS" if okA else "FAIL"))
print("(B) M homomorphic to adjoint(M) via order-2 T (rightremainder = 0, exact)     : %s"
      % ("PASS" if okB else "FAIL"))
print("(C) t=0 Jordan: single 3x3 block, maximal log power n = 2 (even)              : %s"
      % ("PASS" if okC else "FAIL"))
print("(D) det character -1 at {1/9,1/5,1/4,inf} => G = O(3,C); twist u^2=(1-4t)(1-5t)(1-9t): %s"
      % ("PASS" if okD else "FAIL"))
print("---------------------------------------------------------------")
if allpass:
    print("=> The solution space of the IRREDUCIBLE order-3 operator M carries a monodromy-")
    print("   invariant nondegenerate SYMMETRIC bilinear form; the differential Galois group is")
    print("   G = O(3,C), with identity component G^0 = SO(3,C) ~= PSL(2,C).  Hence M IS")
    print("   projectively equivalent to the symmetric square of a second-order operator")
    print("   (defined over Q(t,u), u^2 = (1-4t)(1-5t)(1-9t)), even though M is NOT a LITERAL")
    print("   symmetric square.  The paper's 'excludes every gauge/pullback-equivalent of a")
    print("   symmetric square' and 'SL3-type Galois group / no elliptic closed form' claims")
    print("   are REFUTED: an elliptic/2F1 closed form is EXPECTED, via a 2nd-order operator on")
    print("   that genus-one curve.  Non-Liouvillian (no algebraic/elementary form) still holds,")
    print("   now because G^0 = SO(3,C) is simple hence non-solvable.")
    print("   [Facts (i)-(iii) due to J.-M. Maillard, private communication, July 2026;")
    print("    independently re-verified here in exact arithmetic.]")
else:
    print("=> CERTIFICATE INCOMPLETE -- see FAIL markers above.")
    sys.exit(1)
