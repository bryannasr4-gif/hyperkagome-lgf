"""certify_pullback.py -- the 2F1 pullback of V2 in radicals: exact certificate.

Statement certified (data in pullback_data.json; the paper's variable t is called x below):

  (1) PROOF that the modular function H = 1728/j, as a function on X(Gamma_0(30)+) with the
      coordinate t of Eqs. (10)-(11), satisfies the bidegree-(8,72) polynomial P_H(H,t) = 0
      exactly.  The verification runs through q^1200; it is a proof, not a series check,
      because of the a-priori divisor bound: H lies in the function field of X_0(30), which
      has degree |W| = 8 over Q(t) (W = the Atkin-Lehner group of level 30; covolume ratio
      24pi/3pi = 8), and the pole divisor of H on X_0(30) has degree 72: j has degree
      [PSL_2(Z):Gamma_0(30)] = 72, and nu_3(30) = 0 means no point above j = 0 retains an
      order-three stabiliser, so all 72/3 = 24 of them are ramified of index 3 and 1728/j
      has 24 triple poles.  (nu_3 = 0 does NOT mean j omits the value 0: the map is
      surjective; it means none of its 24 preimages is an elliptic point.)  Hence a
      polynomial G(H,t) of bidegree at most (8,72) not vanishing identically on the curve
      of H has, as a function on that curve, polar degree at most 8*72 + 72*8 = 1152, and
      since a nonzero function on a compact curve has as many zeros as poles it cannot
      vanish at the cusp to order > 1152.  Vanishing through order 1200 forces
      P_H(H,t) == 0.  (q is a uniformiser there: Gamma_0(30) has width one at infinity,
      and ord_q = ord_t because t(q) = q + O(q^2).)
  (2) PROOF that the explicit radical expression H_r (M. van Hoeij, private communication,
      July 2026) is an exact root of P_H: a polynomial identity in the multiquadratic
      algebra A = Q(x)[r1,r2,r3]/(r1^2-(1-x)(1-9x), r2^2-(1-x)(1-5x), r3^2-(1-4x)),
      computed in exact integer arithmetic -- no series involved.
  (3) PROOF that A is a field of degree 8 over Q(x): all seven subset products of the
      radicands are non-squares in Q(x).
  (4) PROOF that the eight sign-conjugates of H_r are pairwise distinct and that the
      analytic branch (signs -1,-1,-1) has the same t-expansion as the modular H.  With
      (1)-(3) this identifies the modular H with H_r exactly and makes P_H irreducible over
      Q(t) (degree 8, annihilating an element of degree 8), hence THE primitive minimal
      polynomial (content of its 657 integer coefficients = 1).  Corollaries: the splitting
      field of P_H is Q(t)(r1,r2,r3) = Q(X_0(30)); Gal(P_H) = (Z/2)^3 = W; the pullback is
      solvable in radicals; and r1 r2 r3 = (1-t) v, v^2 = (1-4t)(1-5t)(1-9t) the twist.
  (5) For the base 2F1([1/8,3/8];[1];.): the bidegree-(4,24) polynomial OtherPH
      annihilates the explicit two-radical expression
      OtherH exactly (radicands (1-x)(1-5x) and (1-x)(1-4x)(1-9x) = r1^2 r3^2: an index-two
      subfield); its four conjugates are distinct (irreducible, minimal); and its analytic
      branch H = 256 t - 25600 t^2 + ... satisfies the Schwarz pullback identity
          Q_E(H) H'^2 + (1/2){H,t} == Q_V(t)
      against the CERTIFIED projective invariant Q_V of V2 (V2_data.json) through t^96
      (guess-and-verify standard; margin stated).

Every primitive is validated on controls of known structure first, and every stage carries
a negative control shown to FAIL.  Writes CERTIFICATE_pullback.txt itself.

    python numerics/certify_pullback.py
"""
import json
import math
import os
import sys
from fractions import Fraction as F

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
LINES = []
FAIL = []


def out(s=""):
    print(s)
    LINES.append(s)


def chk(name, cond, detail=""):
    tag = "PASS" if cond else "FAIL"
    if not cond:
        FAIL.append(name)
    out("  [%s] %s%s" % (tag, name, ("  " + detail) if detail else ""))


# ======================================================================
# integer power series via Kronecker substitution (single big-int multiply)
# ======================================================================
def _encode(a, W):
    n = len(a)
    if n == 1:
        return a[0]
    h = n // 2
    return _encode(a[:h], W) + (_encode(a[h:], W) << (W * h))


def _decode(c, W, n, truncated):
    """balanced-digit decode; digits must satisfy |d| < 2^(W-2)."""
    half = 1 << (W - 1)
    mask = (1 << W) - 1
    o = []
    for _ in range(n):
        r = c & mask
        if r >= half:
            r -= (1 << W)
            c += (1 << W)
        o.append(r)
        c >>= W
    if not truncated:
        assert c == 0, "decode overflow"
    lim = 1 << (W - 2)
    assert all(abs(v) < lim for v in o), "digit width exceeded"
    return o


def kmul(a, b, n, W):
    """(a*b) truncated to n terms; exact for |result digits| < 2^(W-2)."""
    c = _encode(a, W) * _encode(b, W)
    full = len(a) + len(b) - 1
    return _decode(c, W, min(n, full), n < full) + [0] * max(0, n - full)


def ser_inv_int(a, n):
    """inverse of an integer series with a[0] == 1 (integer result)."""
    assert a[0] == 1
    o = [0] * n
    o[0] = 1
    for m in range(1, n):
        s = 0
        for k in range(1, min(m, len(a) - 1) + 1):
            ak = a[k]
            if ak:
                s += ak * o[m - k]
        o[m] = -s
    return o


def ser_div_int(num, den, n):
    """num/den with den[0] == 1, integer quotient, via the recurrence (keeps
    intermediates at the size of the QUOTIENT, avoiding huge inverse series)."""
    assert den[0] == 1
    o = [0] * n
    for m in range(n):
        s = num[m] if m < len(num) else 0
        for k in range(1, min(m, len(den) - 1) + 1):
            dk = den[k]
            if dk:
                s -= dk * o[m - k]
        o[m] = s
    return o


def euler_P(n, step=1):
    """prod_{m>=1} (1 - q^(step*m)), length n (pentagonal, sparse)."""
    o = [0] * n
    o[0] = 1
    k = 1
    while True:
        done = True
        for g in (k * (3 * k - 1) // 2, k * (3 * k + 1) // 2):
            e = g * step
            if e < n:
                o[e] += (-1) ** k
                done = False
        if done:
            break
        k += 1
    return o


# ======================================================================
# multiquadratic algebra: element = {bitmask: integer coeff list in x}
# ======================================================================
def pmul(a, b, W):
    return kmul(a, b, len(a) + len(b) - 1, W)


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def alg_mul(A, B, rads, W):
    o = {}
    for s, ca in A.items():
        for t_, cb in B.items():
            m = s ^ t_
            c = pmul(ca, cb, W)
            common = s & t_
            i = 0
            while common:
                if common & 1:
                    c = pmul(c, rads[i], W)
                common >>= 1
                i += 1
            o[m] = padd(o[m], c) if m in o else c
    return {m: c for m, c in o.items() if any(v != 0 for v in c)}


def alg_add(A, B):
    o = dict(A)
    for m, c in B.items():
        o[m] = padd(o[m], c) if m in o else c
    return {m: c for m, c in o.items() if any(v != 0 for v in c)}


def branch_series(comp, den, rads, signs, n):
    """t-expansion of (sum_masks comp*radicals)/den under given radical signs."""
    g = len(rads)
    rs = []
    for i in range(g):
        rad = [F(c) for c in rads[i]] + [F(0)] * n
        assert rad[0] == 1
        s = [F(0)] * n
        s[0] = F(1)
        for m in range(1, n):
            s[m] = (rad[m] - sum(s[k] * s[m - k] for k in range(1, m))) / 2
        rs.append([F(signs[i]) * c for c in s])

    def fmul(a, b):
        o = [F(0)] * n
        for i in range(n):
            if a[i]:
                for k in range(n - i):
                    if b[k]:
                        o[i + k] += a[i] * b[k]
        return o

    tot = [F(0)] * n
    for mask, c in comp.items():
        term = [F(v) for v in c[:n]] + [F(0)] * max(0, n - len(c))
        for i in range(g):
            if mask >> i & 1:
                term = fmul(term, rs[i])
        tot = [p + q for p, q in zip(tot, term)]
    dfr = [F(c) for c in den[:n]] + [F(0)] * max(0, n - len(den))
    inv = [F(0)] * n
    inv[0] = 1 / dfr[0]
    for m in range(1, n):
        inv[m] = -sum(dfr[k] * inv[m - k] for k in range(1, m + 1)) * inv[0]
    return fmul(tot, inv)


# ======================================================================
out("=" * 78)
out("CERTIFICATE: the 2F1 pullback of V_2 in radicals  (certify_pullback.py)")
out("=" * 78)

D = json.load(open(os.path.join(HERE, "pullback_data.json")))
x = sp.Symbol("x")
RAD_POLYS = [[1, -10, 9], [1, -6, 5], [1, -4]]   # (1-x)(1-9x), (1-x)(1-5x), (1-4x)
W1 = 1408                                        # digit width, Stage 1 (asserted at decode)
W2 = 4096                                        # digit width, Stage 2

out("\n--- Stage 0: validation of primitives on controls -----------------------")
N0 = 20
P1 = euler_P(N0)
P24 = P1
for _ in range(2):
    P24 = kmul(P24, P24, N0, W1)                 # P^4
P8 = kmul(P24, P24, N0, W1)                      # P^8
P16 = kmul(P8, P8, N0, W1)
P24 = kmul(P16, P8, N0, W1)                      # P^24
E4 = [1] + [240 * sum(d ** 3 for d in range(1, n + 1) if n % d == 0) for n in range(1, N0)]
E4c = kmul(kmul(E4, E4, N0, W1), E4, N0, W1)
qj = kmul(E4c, ser_inv_int(P24, N0), N0, W1)     # q*j = E4^3/P^24
chk("j-series control: q*j = 1 + 744 q + 196884 q^2 + ...",
    qj[0] == 1 and qj[1] == 744 and qj[2] == 196884)
chk("E4 head 1 + 240q + 2160q^2", E4[:3] == [1, 240, 2160])
chk("Kronecker-multiply control: (1+q)^2 = 1+2q+q^2 and (1-3q)(1+3q) = 1-9q^2",
    kmul([1, 1], [1, 1], 3, 64) == [1, 2, 1] and kmul([1, -3], [1, 3], 3, 64) == [1, 0, -9])
sq = alg_mul({0b111: [1]}, {0b111: [1]}, RAD_POLYS, W2)
target = [int(c) for c in sp.Poly(sp.expand((1 - x) ** 2 * (1 - 4 * x) * (1 - 5 * x) * (1 - 9 * x)), x).all_coeffs()[::-1]]
got = sq.get(0, [])
chk("algebra control: (r1 r2 r3)^2 == ((1-x) v)^2, v^2 = (1-4x)(1-5x)(1-9x)",
    list(sq.keys()) == [0] and padd(got, [0] * len(target)) == padd(target, [0] * len(got)))
chk("algebra NEGATIVE control: (r1 r2)^2 != ((1-x) v)^2",
    padd(alg_mul({0b011: [1]}, {0b011: [1]}, RAD_POLYS, W2).get(0, []), [0] * len(target))
    != padd(target, [0] * 8))


def is_square_poly(p_expr):
    c, fl = sp.factor_list(sp.expand(p_expr))
    return all(e_ % 2 == 0 for _, e_ in fl) and sp.sqrt(c).is_rational


chk("square-detector control: (2-3x)^2 detected as a square", is_square_poly((2 - 3 * x) ** 2))
chk("square-detector NEGATIVE control: (1-4x) is not a square", not is_square_poly(1 - 4 * x))

# ======================================================================
out("\n--- Stage 1: P_H annihilates the modular H = 1728/j (PROOF) -------------")
N = 1201
# Coefficient growth (radius of convergence in q): t(q) is singular at the W_30
# elliptic point q = e^(-2 pi/sqrt(30)) ~ 0.3176, so |t_n| ~ 10^(0.50 n); H(q) is
# singular at the zero of E4, |q_rho| = e^(-pi sqrt(3)) ~ 0.00433, so |H_n| ~
# 10^(2.36 n) (cross-checked against the observed 2271-digit coefficient at order 962).
W_t = 2560                                       # digit width for t-powers
W_H = 10496                                      # digit width for H-powers
W_F = 13312                                      # digit width for the final products
num = kmul(kmul(euler_P(N, 1), euler_P(N, 6), N, W1), kmul(euler_P(N, 10), euler_P(N, 15), N, W1), N, W1)
den = kmul(kmul(euler_P(N, 2), euler_P(N, 3), N, W1), kmul(euler_P(N, 5), euler_P(N, 30), N, W1), N, W1)
U1 = ser_div_int(num, den, N)
U = kmul(kmul(U1, U1, N, W1), U1, N, W1)         # U = q*u
chk("q*u head [1,-3,3,-1,0,0,0,-3,9]", U[:9] == [1, -3, 3, -1, 0, 0, 0, -3, 9])
U2 = kmul(U, U, N, W1)
denom_t = [U2[i] + (7 * U[i - 1] if i >= 1 else 0) + (1 if i == 2 else 0) for i in range(N)]
tq = ser_div_int([0] + U[:N - 1], denom_t, N)
chk("t(q) head q - 4q^2 + 12q^3 - 34q^4 + 90q^5", tq[:6] == [0, 1, -4, 12, -34, 90])
P1N = euler_P(N)
P4N = kmul(kmul(P1N, P1N, N, W1), kmul(P1N, P1N, N, W1), N, W1)
P8N = kmul(P4N, P4N, N, W1)
P24N = kmul(kmul(P8N, P8N, N, W1), P8N, N, W1)
E4N = [1] + [240 * sum(d ** 3 for d in range(1, n + 1) if n % d == 0) for n in range(1, N)]
E4c3 = kmul(kmul(E4N, E4N, N, W1), E4N, N, W1)
Hq = ser_div_int([0] + [1728 * c for c in P24N[:N - 1]], E4c3, N)
chk("H(q) head 1728 q", Hq[0] == 0 and Hq[1] == 1728)
out("      series built; evaluating P_H(H(q), t(q)) over Z ...")
dH, dx = D["PH"]["dH"], D["PH"]["dx"]
A = [int(c) for c in D["PH"]["coeffs"]]
tpow = [[1] + [0] * (N - 1)]
for k in range(1, dx + 1):
    tpow.append(kmul(tpow[-1], tq, N, W_t))
Hpow = [[1] + [0] * (N - 1)]
for i in range(1, dH + 1):
    Hpow.append(kmul(Hpow[-1], Hq, N, W_H))
total = [0] * N
for i in range(dH + 1):
    row = [0] * N
    for k in range(dx + 1):
        c = A[i * (dx + 1) + k]
        if c:
            tk = tpow[k]
            for m in range(N):
                if tk[m]:
                    row[m] += c * tk[m]
    prod = kmul(row, Hpow[i], N, W_F)
    total = [p + q for p, q in zip(total, prod)]
out("      evaluation done")
chk("P_H(H(q), t(q)) == 0 through q^%d  (all %d orders exactly zero over Z)" % (N - 1, N),
    all(c == 0 for c in total))
out("      a-priori bound: deg(H) = 72 (24 triple poles over j=0, since nu_3(30)=0 leaves no")
out("      elliptic point there) and deg(t) = |W| = 8, so a nonzero G of bidegree <= (8,72)")
out("      has polar degree <= 8*72 + 72*8 = 1152, hence vanishing order <= 1152 at the cusp")
out("      (q is a uniformiser: width one at infinity); 1200 > 1152, so P_H(H,t) == 0 is PROVEN.")
pert = [p + q for p, q in zip(total, tq)]        # simulate A[(i,k)=(0,1)] += 1
chk("NEGATIVE control: corrupting coefficient (i,k)=(0,1) by +1 is detected at q^1",
    pert[1] == 1 and any(c != 0 for c in pert))

# ======================================================================
out("\n--- Stage 2: the radical expression is an EXACT root of P_H --------------")
HR = D["H_radical"]
comp = {}
for mstr, cl in HR["components"].items():
    mask = (1 if mstr[0] == "1" else 0) | (2 if mstr[1] == "1" else 0) | (4 if mstr[2] == "1" else 0)
    comp[mask] = [int(c) for c in cl]
Dden = [int(c) for c in HR["den"]]
Epow = [{0: [1]}]
for i in range(1, dH + 1):
    Epow.append(alg_mul(Epow[-1], comp, RAD_POLYS, W2))
Dpow = [[1]]
for i in range(1, dH + 1):
    Dpow.append(pmul(Dpow[-1], Dden, W2))
TOT = {}
for i in range(dH + 1):
    Ai = [A[i * (dx + 1) + k] for k in range(dx + 1)]
    scal = pmul(Ai, Dpow[dH - i], W2)
    term = {m: pmul(c, scal, W2) for m, c in Epow[i].items()}
    TOT = alg_add(TOT, term)
chk("sum_i A_i(x) Num^i D^(8-i) == 0 in Q(x)[r1,r2,r3]  (all 8 components identically 0)",
    len(TOT) == 0, "exact integer arithmetic, no series")
comp_bad = {m: list(c) for m, c in comp.items()}
comp_bad[0] = padd(comp_bad.get(0, [0]), [1])
Eb = [{0: [1]}]
for i in range(1, dH + 1):
    Eb.append(alg_mul(Eb[-1], comp_bad, RAD_POLYS, W2))
TOTb = {}
for i in range(dH + 1):
    Ai = [A[i * (dx + 1) + k] for k in range(dx + 1)]
    scal = pmul(Ai, Dpow[dH - i], W2)
    TOTb = alg_add(TOTb, {m: pmul(c, scal, W2) for m, c in Eb[i].items()})
chk("NEGATIVE control: perturbing the radical expression by +1 gives a nonzero result",
    len(TOTb) > 0, "exact integer arithmetic, no series")

# ======================================================================
out("\n--- Stage 3: the multiquadratic algebra is a degree-8 field --------------")
rad_exprs = [(1 - x) * (1 - 9 * x), (1 - x) * (1 - 5 * x), (1 - 4 * x)]
names = ["r1^2", "r2^2", "r3^2"]
allok = True
for mask in range(1, 8):
    p = sp.Integer(1)
    nm = []
    for i in range(3):
        if mask >> i & 1:
            p *= rad_exprs[i]
            nm.append(names[i])
    okk = not is_square_poly(p)
    allok = allok and okk
    out("      %-18s non-square: %s" % ("*".join(nm), okk))
chk("all seven subset products of the radicands are non-squares => [A:Q(x)] = 8", allok)

# ======================================================================
out("\n--- Stage 4: conjugates, branch identification, minimality ---------------")
NB = 20      # the slowest-separating conjugate pair first differs at order 15
branches = {}
for sgn in range(8):
    signs = [1 if sgn >> i & 1 else -1 for i in range(3)]
    branches[tuple(signs)] = branch_series(comp, Dden, RAD_POLYS, signs, NB)
keys = list(branches)
chk("the 8 sign-conjugates of the radical expression are pairwise distinct",
    all(branches[keys[i]] != branches[keys[j]]
        for i in range(8) for j in range(i + 1, 8)))
# modular H as a t-series: revert t(q), compose H(q) with q(t)
rev = [F(0)] * NB
rev[1] = F(1)


def comp_small(series, inner, n):
    o = [F(0)] * n
    pw = [F(1)] + [F(0)] * (n - 1)
    o[0] = F(series[0])
    for k in range(1, n):
        npw = [F(0)] * n
        for a_ in range(n):
            if pw[a_]:
                for b_ in range(n - a_):
                    if inner[b_]:
                        npw[a_ + b_] += pw[a_] * inner[b_]
        pw = npw
        if series[k]:
            o = [p + F(series[k]) * q for p, q in zip(o, pw)]
    return o


for _ in range(NB + 2):
    errv = comp_small(tq[:NB], rev, NB)
    errv[1] -= 1
    rev = [rev[i] - errv[i] for i in range(NB)]
chk("series reversion control: t(q(t)) == t through t^%d" % (NB - 1),
    comp_small(tq[:NB], rev, NB) == [F(0), F(1)] + [F(0)] * (NB - 2))
Hx = comp_small(Hq[:NB], rev, NB)
chk("modular H(t) head 1728t - 1278720t^2 + 606044160t^3",
    Hx[1] == 1728 and Hx[2] == -1278720 and Hx[3] == 606044160)
anal = branches[(-1, -1, -1)]
chk("the analytic branch (-1,-1,-1) equals the modular H(t) through t^%d" % (NB - 1),
    all(anal[i] == Hx[i] for i in range(NB)))
chk("NEGATIVE control: every other sign branch differs from H(t) within %d orders" % NB,
    all(any(branches[s][i] != Hx[i] for i in range(NB))
        for s in branches if s != (-1, -1, -1)))
g = 0
for c in A:
    g = math.gcd(g, abs(c))
    if g == 1:
        break
chk("content of the 657 integer coefficients of P_H is 1 (primitive)", g == 1)
chk("deg_H P_H = 8 (leading block nonzero)", any(A[8 * (dx + 1) + k] for k in range(dx + 1)))
chk("deg_t P_H = 72 (some t^72 coefficient nonzero)", any(A[i * (dx + 1) + 72] for i in range(dH + 1)))
out("      Conclusion: P_H (degree 8) annihilates an element of exact degree 8, so it is")
out("      irreducible over Q(t) and IS the primitive minimal polynomial, bidegree (8,72).")
out("      Its splitting field is Q(t)(r1,r2,r3) = Q(X_0(30)); Gal = (Z/2)^3 = the")
out("      Atkin-Lehner group of level 30; the pullback is solvable in radicals.")

# ======================================================================
out("\n--- Stage 5: the (4,24) equation for the 2F1([1/8,3/8];[1];.) base -------")
OP = D["OtherPH"]
oco = {int(i): [int(c) for c in cl] for i, cl in OP["coeffs_by_H_power"].items()}
RAD2 = [[1, -6, 5], [1, -14, 49, -36]]           # (1-x)(1-5x), (1-x)(1-4x)(1-9x)
OH = D["OtherH"]
ocomp = {}
for mstr, cl in OH["components"].items():
    mask = (1 if mstr[0] == "1" else 0) | (2 if mstr[1] == "1" else 0)
    ocomp[mask] = [int(c) for c in cl]
oden = [int(c) for c in OH["den"]]
E2p = [{0: [1]}]
for i in range(1, 5):
    E2p.append(alg_mul(E2p[-1], ocomp, RAD2, W2))
D2p = [[1]]
for i in range(1, 5):
    D2p.append(pmul(D2p[-1], oden, W2))
TOT2 = {}
for i in range(5):
    scal = pmul(oco.get(i, [0]), D2p[4 - i], W2)
    TOT2 = alg_add(TOT2, {m: pmul(c, scal, W2) for m, c in E2p[i].items()})
chk("OtherPH(OtherH) == 0 exactly in Q(x)[r2,s2]  (all 4 components identically 0)",
    len(TOT2) == 0)
ob = {}
for sgn in range(4):
    signs = [1 if sgn >> i & 1 else -1 for i in range(2)]
    ob[tuple(signs)] = branch_series(ocomp, oden, RAD2, signs, 8)
k2 = list(ob)
chk("its 4 sign-conjugates are pairwise distinct => irreducible => minimal, bidegree (4,24)",
    all(ob[k2[i]] != ob[k2[j]] for i in range(4) for j in range(i + 1, 4)))
chk("s2^2 == r1^2 * r3^2 / (1-x)^0  i.e. (1-x)(1-4x)(1-9x) == [(1-x)(1-9x)]*(1-4x)",
    sp.expand((1 - x) * (1 - 4 * x) * (1 - 9 * x) - rad_exprs[0] * rad_exprs[2]) == 0,
    "so s2 = +-r1 r3 and Q(x)(r2, s2) is an index-two subfield of Q(x)(r1,r2,r3)")
# Schwarz identity against certified Q_V
NS = 101
V2 = json.load(open(os.path.join(HERE, "V2_data.json")))
tsym = sp.Symbol("t")
QV = sp.cancel(sp.sympify(V2["Q"]["num"]) / sp.sympify(V2["Q"]["den"]))
z = sp.Symbol("z")
P_base = sp.cancel((sp.Rational(3, 2) * z - 1) / (z * (z - 1)))
Q_base = sp.cancel(sp.Rational(3, 64) / (z * (z - 1)))
Q_E = sp.cancel(Q_base - sp.diff(P_base, z) / 2 - P_base ** 2 / 4)
chk("Q_E of the printed base operator = (15z^2-19z+16)/(64 z^2 (z-1)^2)",
    sp.cancel(Q_E - (15 * z ** 2 - 19 * z + 16) / (64 * z ** 2 * (z - 1) ** 2)) == 0)


def qmulF(a, b, n):
    o = [F(0)] * n
    for i in range(n):
        if a[i]:
            for k in range(n - i):
                if b[k]:
                    o[i + k] += a[i] * b[k]
    return o


def qinvF(a, n):
    o = [F(0)] * n
    o[0] = 1 / a[0]
    for m in range(1, n):
        o[m] = -sum(a[k] * o[m - k] for k in range(1, m + 1)) * o[0]
    return o


def qratF(expr, var, n):
    nu, de = sp.fraction(sp.cancel(sp.together(expr)))
    ns = [F(int(c)) for c in sp.Poly(nu, var).all_coeffs()[::-1]] + [F(0)] * n
    ds = [F(int(c)) for c in sp.Poly(de, var).all_coeffs()[::-1]] + [F(0)] * n
    return qmulF(ns[:n], qinvF(ds[:n], n), n)


CO5 = []
for k in range(5):
    cs = oco.get(k, [0])
    shift = 4 - k
    o5 = [F(0)] * NS
    for e_, c in enumerate(cs):
        if e_ - shift >= 0:
            if e_ - shift < NS:
                o5[e_ - shift] = F(c)
        else:
            assert c == 0
    CO5.append(o5)


def PE5(kk, dv=False):
    pw = [[F(0)] * NS for _ in range(5)]
    pw[0][0] = F(1)
    for i2 in range(1, 5):
        pw[i2] = qmulF(pw[i2 - 1], kk, NS)
    acc = [F(0)] * NS
    for i2 in range(5):
        if dv:
            if i2 >= 1:
                acc = [p + F(i2) * q for p, q in zip(acc, qmulF(CO5[i2], pw[i2 - 1], NS))]
        else:
            acc = [p + q for p, q in zip(acc, qmulF(CO5[i2], pw[i2], NS))]
    return acc


kk = [F(0)] * NS
kk[0] = F(256)
for _ in range(9):
    kk = [a - b for a, b in zip(kk, qmulF(PE5(kk), qinvF(PE5(kk, dv=True), NS), NS))]
chk("Newton branch of OtherPH exact to t^%d" % (NS - 1), all(c == 0 for c in PE5(kk)))
h5 = [F(0)] + kk[:NS - 1]
chk("branch head 256t - 25600t^2 + 1441792t^3",
    h5[1] == 256 and h5[2] == -25600 and h5[3] == 1441792)


def dser5(a):
    return [a[i + 1] * (i + 1) for i in range(NS - 1)] + [F(0)]


def schwarz_residual(h):
    h1 = dser5(h)
    h2 = dser5(h1)
    h3 = dser5(h2)
    ih1 = qinvF(h1, NS)
    schw = [p - F(3, 2) * q for p, q in
            zip(qmulF(h3, ih1, NS), qmulF(qmulF(h2, ih1, NS), qmulF(h2, ih1, NS), NS), )]
    qe_reg = qratF(sp.cancel(Q_E * z ** 2), z, NS)
    qe = [F(0)] * NS
    pwc = [F(1)] + [F(0)] * (NS - 1)
    for k in range(NS):
        if qe_reg[k]:
            qe = [p + qe_reg[k] * q for p, q in zip(qe, pwc)]
        if k < NS - 1:
            pwc = qmulF(pwc, h, NS)
    hoverT = [h[i + 1] for i in range(NS - 1)] + [F(0)]
    fac = qmulF(h1, qinvF(hoverT, NS), NS)
    lhs1 = qmulF(qe, qmulF(fac, fac, NS), NS)
    lhs2 = [F(0), F(0)] + [F(1, 2) * schw[i] for i in range(NS - 2)]
    rhs = qratF(sp.cancel(QV * tsym ** 2), tsym, NS)
    return [lhs1[i] + lhs2[i] - rhs[i] for i in range(NS - 4)]


res5 = schwarz_residual(h5)
chk("Schwarz identity Q_E(H)H'^2 + (1/2){H,t} == Q_V through t^%d" % (NS - 5),
    all(c == 0 for c in res5),
    "against certified V2_data.json (guess-and-verify standard, %d orders)" % (NS - 4))
hbad = list(h5)
hbad[6] += F(1)
chk("NEGATIVE control: perturbing the branch at t^6 breaks the Schwarz identity",
    any(c != 0 for c in schwarz_residual(hbad)))

# ======================================================================
out("\n" + "=" * 78)
ok = len(FAIL) == 0
out("RESULT: %s  (%d checks failed)"
    % ("ALL PASS" if ok else "FAILURE", len(FAIL)))
if not ok:
    for f in FAIL:
        out("  FAILED: " + f)
out("The 2F1([1/12,5/12];[1];.) pullback H = 1728/j of V_2 has primitive minimal")
out("polynomial P_H of bidegree (8,72) over Q(t); its splitting field is")
out("Q(t)(sqrt((1-t)(1-9t)), sqrt((1-t)(1-5t)), sqrt(1-4t)) = Q(X_0(30)), with Galois")
out("group (Z/2)^3, the Atkin-Lehner group of level 30; the product of the three")
out("radicands is ((1-t)v)^2, v the determinant-character twist.  The pullback is")
out("solvable in radicals (explicit expression: M. van Hoeij, private communication,")
out("July 2026).  For the 2F1([1/8,3/8];[1];.) base the minimal polynomial has")
out("bidegree (4,24) on an index-two subfield.")
out("=" * 78)
open(os.path.join(HERE, "CERTIFICATE_pullback.txt"), "w").write("\n".join(LINES) + "\n")
print("wrote", os.path.join(HERE, "CERTIFICATE_pullback.txt"))
sys.exit(0 if ok else 1)
