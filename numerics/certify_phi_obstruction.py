"""certify_phi_obstruction.py -- Phi itself is NOT weight-two quasimodular (theorem),
and the exact identity that isolates the obstruction.

Setting (established elsewhere in this repository):
  t = the degree-one coordinate on X(Gamma_0(30)+) of certify_modular.py,
  W = q dt/dq,  ' = d/dt,  v = sqrt((1-4t)(1-5t)(1-9t)), v(0)=+1,
  y0 = Phi'/2 = [rho0 W + rho1 W']/v          (certify_y0.py),
  Q_V from V2_data.json;  Ntilde = D^3 + 4 Q_V D + 2 Q_V' annihilates W,
  Gal(w'' + Q_V w = 0) = SL(2,C)   (certified in certify_y0_lemma.py).
  Two consequences of that lemma are used here; both are short exact algebra:
    (i)  Ntilde is IRREDUCIBLE (its solution space is the symmetric square of
         the standard SL_2 representation, which is irreducible), so no
         operator of order <= 2 annihilates W: W, W', W'' are linearly
         independent over Q(t).
    (ii) no nonzero element  xi = S0 W + S1 W' + S2 W''  (S_i in Q(t)) is
         algebraic over Q(t): SL(2,C) is connected, so the algebraic closure
         of C(t) in the Picard-Vessiot field is C(t) itself and xi would be
         RATIONAL; then P = S0 + S1 D + S2 D^2 sends W to xi in C(t), the
         Galois orbit of W spans Sol(Ntilde) (irreducibility), so P is
         constant on a spanning orbit and kills a >= 2-dimensional subspace
         of Sol(Ntilde); GCRD(P, Ntilde) would be a proper right factor of
         the irreducible Ntilde unless P = 0, whence xi = 0.

THIS CERTIFICATE PROVES, in that setting:

  (A) IDENTITY   Phi = 2 rho1 W/v + 2 Int Delta(t) W/v dt + 2/15,
                 Delta := rho0 - (rho1' - rho1 v'/v) = (13t-1) / (30 t (t-1)^2) != 0.

  (B) REDUCTION  if Phi = C(t) + [S0 W + S1 W' + S2 W'']/v with C, S_i in Q(t),
                 then C' = 0, S1 = -d_v(S2), S0 = 2 rho1 + d_v^2(S2) + 4 Q_V S2, and
                 Ntilde_v(S2) = 2 Delta, where d_v = d/dt - v'/v and
                 Ntilde_v = v Ntilde v^{-1} = d_v^3 + 4 Q_V d_v + 2 Q_V'.
                 (verified symbolically for a GENERIC S2; C'=0 uses (ii):
                 a nonzero rational C' would equal an element of the module.)

  (C) THEOREM    Ntilde_v(S2) = 2 Delta has NO rational solution:
                 the Abramov denominator bound at each singular point gives
                 m_max = 0 (no poles possible), the degree bound at infinity gives
                 deg <= 3, and the resulting 4-dimensional linear system is
                 inconsistent, exactly over Q.
    ==>  Phi is NOT in  Q(t) + v^{-1} span_{Q(t)}{W, W', W''};
         since that module contains every chi-twisted meromorphic quasimodular
         form of weight two and depth <= 2, Phi -- unlike its derivative -- is
         not weight-two quasimodular.  The one remaining transcendental in Phi
         is the Eichler-type integral of the weight-four form Delta W^2 / v.

Physically:  (E-1) G_disp(E) = Phi(t), t = (E-1)^{-2}, so the modular closed form
of y0 determines  d/dE [ (E-1) G_disp(E) ] = -4 y0(t) / (E-1)^3  exactly, while
(E-1) G_disp itself is provably outside the weight-two quasimodular module.

Method notes.  Every primitive is validated on controls with known outcomes
before the real run, and every check is paired with a falsification (a perturbed
input that must FAIL).  Series arithmetic is exact rational; the rational-
solution step is exact linear algebra whose sampling count exceeds the a-priori
degree of the polynomial identity, so inconsistency is a proof, not numerics.

Run:  python numerics/certify_phi_obstruction.py     (writes CERTIFICATE_phi_obstruction.txt)
"""
import json
import os
import sys
from fractions import Fraction as F

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
t = sp.symbols("t")
al = sp.symbols("alpha")

ORD = 200          # orders checked for every q-series identity
NQ = ORD + 10

LINES = []


def emit(s=""):
    print(s)
    LINES.append(s)


ok_all = True


def check(name, cond, detail=""):
    global ok_all
    ok_all = ok_all and bool(cond)
    emit("  [%s] %s%s" % ("PASS" if cond else "FAIL", name,
                          ("   " + detail) if detail else ""))
    return cond


# =====================================================================
# series toolkit (exact Fractions, ascending in q) -- as in certify_y0.py
# =====================================================================
def smul(a, b, n=NQ):
    r = [F(0)] * n
    for i, ai in enumerate(a[:n]):
        if ai:
            for j, bj in enumerate(b[:n - i]):
                if bj:
                    r[i + j] += ai * bj
    return r


def sinv(a, n=NQ):
    inv = F(1) / a[0]
    r = [F(0)] * n
    r[0] = inv
    for m in range(1, n):
        acc = F(0)
        for k in range(1, m + 1):
            ak = a[k] if k < len(a) else F(0)
            if ak:
                acc += ak * r[m - k]
        r[m] = -acc * inv
    return r


def sdiv(a, b, n=NQ):
    return smul(a, sinv(b, n), n)


def ssqrt(a, n=NQ):
    assert a[0] == 1
    s = [F(0)] * n
    s[0] = F(1)
    for m in range(1, n):
        acc = a[m] if m < len(a) else F(0)
        for k in range(1, m):
            acc -= s[k] * s[m - k]
        s[m] = acc / 2
    return s


def theta(s, n=NQ):
    return [F(k) * s[k] for k in range(n)]


def scomp(f, g, n=NQ):
    assert g[0] == 0
    acc = [F(0)] * n
    for k in range(min(len(f), n) - 1, -1, -1):
        acc = smul(acc, g, n)
        acc[0] += f[k]
    return acc


def polyser(ci, ts, n=NQ):
    out = [F(0)] * n
    tp = [F(1)] + [F(0)] * (n - 1)
    for k, c in enumerate(ci):
        if k > 0:
            tp = smul(tp, ts, n)
        if c:
            for i in range(n):
                if tp[i]:
                    out[i] += c * tp[i]
    return out


def eta_ipart(d, n=NQ):
    s = [F(0)] * n
    s[0] = F(1)
    k = 1
    while d * k < n:
        new = s[:]
        for i in range(n - 1, d * k - 1, -1):
            new[i] -= s[i - d * k]
        s = new
        k += 1
    return s


emit("certify_phi_obstruction.py -- Phi (unlike Phi') is not weight-two quasimodular")
emit("")
emit("Dependencies: certify_y0.py (the rho-gauge), certify_y0_lemma.py")
emit("(Gal(w''+Q_V w=0) = SL(2,C); irreducibility of Ntilde and the no-algebraic-")
emit("element property of span{W,W',W''} follow -- see the docstring),")
emit("certify_modular.py (t, W), V2_data.json (Q_V), M_coeffs.json + nu.json.")
emit("")

# =====================================================================
# Stage 0: validate the series primitives on known values
# =====================================================================
emit("Stage 0: validate series primitives")
sq = ssqrt([F(1), F(2), F(1)] + [F(0)] * (NQ - 3))
check("sqrt(1+2q+q^2) == 1+q", sq[:3] == [F(1), F(1), F(0)] and
      all(z == 0 for z in sq[3:40]))
inv = sinv([F(1), F(-1)] + [F(0)] * (NQ - 2))
check("1/(1-q) == sum q^k", all(inv[k] == 1 for k in range(40)))
bad = ssqrt([F(1), F(2), F(2)] + [F(0)] * (NQ - 3))
check("NEGATIVE CONTROL: sqrt(1+2q+2q^2) is NOT 1+q", any(z != 0 for z in bad[2:40]))
emit("")

# =====================================================================
# Stage 1: the q-pipeline and the certified gauge identity
# =====================================================================
emit("Stage 1: eta pipeline, y0 from the moments, and the rho-gauge anchor")
Pn = [F(1)] + [F(0)] * (NQ - 1)
for d in (1, 6, 10, 15):
    Pn = smul(Pn, eta_ipart(d))
Pd = [F(1)] + [F(0)] * (NQ - 1)
for d in (2, 3, 5, 30):
    Pd = smul(Pd, eta_ipart(d))
ratio = sdiv(Pn, Pd)
U = smul(smul(ratio, ratio), ratio)                        # q*u
check("q*u head = [1,-3,3,-1,0,0,0,-3,9]",
      [str(x) for x in U[:9]] == ['1', '-3', '3', '-1', '0', '0', '0', '-3', '9'])

qU = [F(0)] + U[:NQ - 1]
q2 = [F(0)] * NQ
q2[2] = F(1)
U2 = smul(U, U)
tser = sdiv(qU, [U2[i] + 7 * qU[i] + q2[i] for i in range(NQ)])   # t = u/(u^2+7u+1)
check("t(q) head = q -4q^2 +12q^3 -34q^4 +90q^5",
      [str(tser[i]) for i in range(1, 6)] == ['1', '-4', '12', '-34', '90'])

W = theta(tser)
prod_poly = [F(1)]
for c1 in (-4, -5, -9):
    prod_poly = [((prod_poly[i] if i < len(prod_poly) else F(0)) +
                  (c1 * prod_poly[i - 1] if i >= 1 else F(0)))
                 for i in range(len(prod_poly) + 1)]
prod_q = scomp(prod_poly, tser)
V = ssqrt(prod_q)
check("v(0)=+1 and v^2 == (1-4t)(1-5t)(1-9t) o t(q) to %d orders" % ORD,
      V[0] == 1 and all(smul(V, V)[i] == prod_q[i] for i in range(ORD)))

Mj = json.load(open(os.path.join(HERE, "M_coeffs.json")))
assert Mj["order"] == 3
Mc = [[F(int(x)) for x in lst] for lst in Mj["coeffs"]]
nu = [F(p, qd) for p, qd in json.load(open(os.path.join(HERE, "nu.json")))]
y0_head = [F(k + 1) * nu[k + 1] / 2 for k in range(len(nu) - 1)]
check("y0 head = [1,10,87,724]", [str(x) for x in y0_head[:4]] == ['1', '10', '87', '724'])


def y0_extend(known, N, c=Mc):
    a = list(known) + [None] * (N - len(known))
    def ffl(n, i):
        r = F(1)
        for k in range(i):
            r *= (n - k)
        return r
    ci = [dict(enumerate(cc)) for cc in c]
    for n0 in range(len(known), N):
        m = n0 - 1
        cu = F(0)
        rhs = F(0)
        for i in range(4):
            for j, cij in ci[i].items():
                if cij == 0:
                    continue
                n = m + i - j
                if n < 0:
                    continue
                val = cij * ffl(n, i)
                if n == n0:
                    cu += val
                elif n < n0:
                    rhs += val * a[n]
        a[n0] = -rhs / cu
    return a


y0_t = y0_extend(y0_head, NQ)
Phi_t = [nu[0]] + [F(2) * y0_t[m - 1] / m for m in range(1, NQ)]
check("Phi head equals the moment sequence nu",
      all(Phi_t[k] == nu[k] for k in range(min(20, len(nu)))))

# the certified gauge (literals identical to certify_y0.py)
B_INT = [0, 0, 30, -600, 4140, -12000, 13830, -5400]
A0_INT = [-8, 231, -2130, 8000, -11436, 2445, 4050]
A1_INT = [0, 8, -169, 1260, -3986, 4432, 1155, -2700]
y0_q = scomp([F(x) for x in y0_t], tser)
thW = theta(W)
lhs = smul(smul(polyser([F(c) for c in B_INT], tser), y0_q), smul(V, W))
rhs = [smul(polyser([F(c) for c in A0_INT], tser), smul(W, W))[i] +
       smul(polyser([F(c) for c in A1_INT], tser), thW)[i] for i in range(NQ)]
check("gauge anchor: B (y0 o t) v W == A0 W^2 + A1 (q dW/dq) to %d orders" % ORD,
      all(lhs[i] == rhs[i] for i in range(ORD)))
emit("")

# =====================================================================
# Stage 2: Delta, and the reduction -- symbolic, generic S2
# =====================================================================
emit("Stage 2: Delta exactly, and the reduction to Ntilde_v(S2) = 2 Delta")
Bp = sum(sp.Integer(c) * t**k for k, c in enumerate(B_INT))
rho0 = sp.cancel(sum(sp.Integer(c) * t**k for k, c in enumerate(A0_INT)) / Bp)
rho1 = sp.cancel(sum(sp.Integer(c) * t**k for k, c in enumerate(A1_INT)) / Bp)
vv = (1 - 4 * t) * (1 - 5 * t) * (1 - 9 * t)
lam = sp.cancel(sp.diff(vv, t) / (2 * vv))
Delta = sp.cancel(rho0 - (sp.diff(rho1, t) - rho1 * lam))
check("Delta == (13t-1)/(30 t (t-1)^2)",
      sp.simplify(Delta - (13 * t - 1) / (30 * t * (t - 1) ** 2)) == 0)
check("Delta != 0  (the S2 = 0 case is impossible at every degree)", Delta != 0)

V2 = json.load(open(os.path.join(HERE, "V2_data.json")))
QV = sp.cancel(sp.sympify(V2["Q"]["num"]) / sp.sympify(V2["Q"]["den"]))

# generic-S2 reduction: represent W, W', W'' as symbols w0,w1,w2 with the
# rewrite d/dt: w0->w1, w1->w2, w2 -> -4 Q_V w1 - 2 Q_V' w0  (Ntilde(W) = 0).
w0, w1, w2 = sp.symbols("w0 w1 w2")
S2f = sp.Function("S2")


def ddt(expr):
    """total t-derivative on Q(t)[w0,w1,w2] under the rewrite rule"""
    e = sp.diff(expr, t)
    e += sp.diff(expr, w0) * w1 + sp.diff(expr, w1) * w2
    e += sp.diff(expr, w2) * (-4 * QV * w1 - 2 * sp.diff(QV, t) * w0)
    return sp.expand(e)


dv = lambda f: sp.diff(f, t) - lam * f
S2 = S2f(t)
S1 = -dv(S2)
S0 = 2 * rho1 + dv(dv(S2)) + 4 * QV * S2
expr = (S0 * w0 + S1 * w1 + S2 * w2)          # this is v * (candidate Phi - C)
# d/dt[ expr / v ] = [ ddt(expr) - lam*expr ] / v   must equal 2 y0 = [2rho0 w0 + 2rho1 w1]/v
resid = sp.expand(ddt(expr) - lam * expr - (2 * rho0 * w0 + 2 * rho1 * w1))
cw0 = sp.simplify(resid.coeff(w0, 1).coeff(w1, 0).coeff(w2, 0))
cw1 = sp.simplify(resid.coeff(w1, 1).coeff(w0, 0).coeff(w2, 0))
cw2 = sp.simplify(resid.coeff(w2, 1).coeff(w0, 0).coeff(w1, 0))
NtvS2 = sp.expand(dv(dv(dv(S2))) + 4 * QV * dv(S2) + 2 * sp.diff(QV, t) * S2)
check("reduction (generic S2): W'' and W' coefficients vanish identically",
      cw1 == 0 and cw2 == 0)
check("reduction (generic S2): W coefficient == Ntilde_v(S2) - 2 Delta",
      sp.simplify(cw0 - (NtvS2 - 2 * Delta)) == 0)
check("NEGATIVE CONTROL: with S1 = +d_v(S2) the W'' coefficient does NOT vanish",
      sp.simplify(sp.expand(ddt(S0 * w0 + dv(S2) * w1 + S2 * w2)
                            - lam * (S0 * w0 + dv(S2) * w1 + S2 * w2)
                            - (2 * rho0 * w0 + 2 * rho1 * w1)).coeff(w2, 1)
                  .coeff(w0, 0).coeff(w1, 0)) != 0)
emit("  (C' = 0 in the reduction: a nonzero rational C' would equal an element of")
emit("   v^{-1} span{W,W',W''}, impossible by consequence (ii) of the SL(2,C) lemma")
emit("   -- see the docstring for the GCRD argument.)")
emit("")

# =====================================================================
# Stage 3: the constructive identity for Phi, with the constant 2/15
# =====================================================================
emit("Stage 3: Phi == 2 rho1 W/v + 2 Int Delta W/v dt + 2/15   to %d orders" % ORD)
Phi_q = scomp(Phi_t, tser)
# 2 rho1 W / v: numerator (2 A1 o t) W, denominator (B o t) v; both start at q^2
numr = smul(polyser([F(2 * c) for c in A1_INT], tser), W)
denr = smul(polyser([F(c) for c in B_INT], tser), V)
assert numr[0] == numr[1] == 0 and denr[0] == denr[1] == 0
term1 = sdiv(numr[2:] + [F(0), F(0)], denr[2:] + [F(0), F(0)])
# Int Delta W/v dt: q A'(q) = Delta(t) W^2 / v; Delta = (13t-1)/(30 t (t-1)^2)
dnum = polyser([F(-1), F(13)], tser)                      # 13t - 1
dden = polyser([F(0), F(30), F(-60), F(30)], tser)        # 30 t (t-1)^2
integ_num = smul(dnum, smul(W, W))                        # order >= 2  (W^2 ~ q^2)
integ_den = smul(dden, V)                                 # order exactly 1 (30t ~ 30q)
assert integ_num[0] == integ_num[1] == 0
assert integ_den[0] == 0 and integ_den[1] != 0
integrand = sdiv(integ_num[1:] + [F(0)], integ_den[1:] + [F(0)])   # = q A'(q), ord >= 1
check("integrand q A'(q) = Delta W^2/v has no constant term (integrable)",
      integrand[0] == 0)
Aq = [F(0)] * NQ
for m in range(1, NQ):
    Aq[m] = integrand[m] / m
resid3 = [Phi_q[i] - term1[i] - 2 * Aq[i] for i in range(NQ)]
check("residual is the CONSTANT 2/15: q^0 coefficient", resid3[0] == F(2, 15),
      "got %s" % resid3[0])
bad3 = next((i for i in range(1, ORD) if resid3[i] != 0), None)
check("residual constant through q^%d (first nonzero = %s)" % (ORD - 1, bad3),
      bad3 is None)
# falsification: Delta with 13 -> 14 must break it
dnum_bad = polyser([F(-1), F(14)], tser)
integ_bad = sdiv(smul(dnum_bad, smul(W, W))[1:] + [F(0)], integ_den[1:] + [F(0)])
Ab = [F(0)] * NQ
for m in range(1, NQ):
    Ab[m] = integ_bad[m] / m
check("NEGATIVE CONTROL: perturbing Delta (13->14) breaks the identity",
      any(Phi_q[i] - term1[i] - 2 * Ab[i] != 0 for i in range(1, 40)))
emit("")

# =====================================================================
# Stage 4: THEOREM -- no rational S2 (Abramov bounds, then exact solve)
# =====================================================================
emit("Stage 4: no rational solution of Ntilde_v(S2) = 2 Delta")


def ffal(n, i):
    r = sp.Integer(1)
    for k in range(i):
        r *= (n - k)
    return r


c3s = sp.Integer(1)
c2s = sp.cancel(-3 * lam)
c1s = sp.cancel(3 * (lam**2 - sp.diff(lam, t)) + 4 * QV)
c0s = sp.cancel(3 * lam * sp.diff(lam, t) - sp.diff(lam, t, 2) - lam**3
                - 4 * QV * lam + 2 * sp.diff(QV, t))
COEF = [c0s, c1s, c2s, c3s]
gs = sp.cancel(2 * Delta)

Acl = sp.Integer(1)
for d in [sp.factor(sp.fraction(sp.cancel(c))[1]) for c in COEF] + [sp.fraction(gs)[1]]:
    Acl = sp.lcm(Acl, d)
CH = []
for c in COEF:
    numd, dend = sp.fraction(sp.together(sp.cancel(Acl * c)))
    assert dend == 1
    CH.append(sp.Poly(sp.expand(numd), t))
GH = sp.Poly(sp.expand(sp.cancel(Acl * gs)), t)
check("denominators cleared; leading coefficient = %s" % sp.factor(Acl), True)


def ordp(poly, p):
    if poly.is_zero:
        return None
    qq = sp.Poly(poly.as_expr().subs(t, t + p), t)
    co = qq.all_coeffs()[::-1]
    for k, c in enumerate(co):
        if c != 0:
            return k, c
    return None


def local_ind(chats, p):
    o = []
    for i, ci in enumerate(chats):
        r = ordp(ci, p)
        o.append((r[0] - i, r[1]) if r else (None, None))
    nu_p = min(vd for vd, _ in o if vd is not None)
    ind = sp.expand(sum(lc * ffal(al, i) for i, (vd, lc) in enumerate(o) if vd == nu_p))
    return nu_p, ind


def pos_roots(ind):
    return [int(r) for r, _ in sp.roots(sp.Poly(ind.subs(al, -al), al)).items()
            if r.is_integer and r > 0]


def bounds(chats, ghat):
    mm = {}
    for p in sp.roots(chats[3]):
        nu_p, ind = local_ind(chats, p)
        w = ordp(ghat, p)
        wp = w[0] if w else 10**9
        mm[p] = max([nu_p - wp] + pos_roots(ind) + [0])
    mu = max(chats[i].degree() - i for i in range(4) if not chats[i].is_zero)
    ind_oo = sp.expand(sum(chats[i].LC() * ffal(al, i) for i in range(4)
                           if not chats[i].is_zero and chats[i].degree() - i == mu))
    iroots = [int(r) for r, _ in sp.roots(sp.Poly(ind_oo, al)).items() if r.is_integer]
    kmax = max([(-10**9 if ghat.is_zero else ghat.degree() - mu)] + iroots)
    return mm, max(kmax + sum(mm.values()), -1), mu


def solve_box(chats, ghat, mm, Nmax):
    Dstar = sp.prod([(t - p) ** m for p, m in mm.items()])
    degD = sum(mm.values())
    if Nmax < 0:
        return None, "degree bound < 0"
    Dmax = max(max((chats[i].degree() if not chats[i].is_zero else 0)
                   + Nmax + i * max(degD - 1, 0) + (3 - i) * degD for i in range(4)),
               (0 if ghat.is_zero else ghat.degree()) + 4 * degD)
    invD = 1 / Dstar
    gders = [invD]
    for _ in range(3):
        gders.append(sp.cancel(sp.diff(gders[-1], t)))
    Hnd = [sp.fraction(sp.together(sp.cancel(
        sum(chats[i].as_expr() * sp.binomial(i, j) * gders[i - j] for i in range(j, 4)))))
        for j in range(4)]
    pts = []
    x = F(2)
    badp = {sp.Rational(p) for p in mm}
    while len(pts) < Dmax + 6:
        if all(sp.Rational(x.numerator, x.denominator) != b for b in badp):
            pts.append(x)
        x += 1

    def ev(nd, x):
        xr = sp.Rational(x.numerator, x.denominator)
        nv = sp.Rational(sp.Poly(nd[0], t).eval(xr))
        dv_ = sp.Rational(sp.Poly(nd[1], t).eval(xr))
        return F(int(nv.p), int(nv.q)) / F(int(dv_.p), int(dv_.q))

    m = []
    for x in pts:
        hv = [ev(nd, x) for nd in Hnd]
        row = []
        for k in range(Nmax + 1):
            s = F(0)
            for j in range(4):
                if k - j >= 0:
                    c = 1
                    for i in range(j):
                        c *= (k - i)
                    s += hv[j] * c * x ** (k - j)
            row.append(s)
        xr = sp.Rational(x.numerator, x.denominator)
        gv = sp.Rational(GH_eval(ghat, xr))
        row.append(F(int(gv.p), int(gv.q)))
        m.append(row)
    nc = Nmax + 1
    r = 0
    piv = []
    for c in range(nc):
        pz = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if pz is None:
            continue
        m[r], m[pz] = m[pz], m[r]
        invv = m[r][c]
        m[r] = [z / invv for z in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                fz = m[i][c]
                m[i] = [m[i][k] - fz * m[r][k] for k in range(nc + 1)]
        piv.append(c)
        r += 1
    if any(all(m[i][c] == 0 for c in range(nc)) and m[i][nc] != 0 for i in range(len(m))):
        return None, "inconsistent (%d samples > identity degree %d)" % (len(pts), Dmax)
    sol = [F(0)] * nc
    for i, c in enumerate(piv):
        sol[c] = m[i][nc]
    N = sum(sp.Rational(sol[k].numerator, sol[k].denominator) * t**k for k in range(nc))
    return sp.cancel(N / Dstar), "candidate found"


def GH_eval(ghat, xr):
    return sp.Integer(0) if ghat.is_zero else ghat.eval(xr)


def applyL(chats, S):
    return sp.cancel(sum(chats[i].as_expr() * sp.diff(S, t, i) for i in range(4)))


# controls first
S_star = (3 * t**2 - 2) / (t**2 * (t - 1) * (4 * t - 1))
gnum, gden = sp.fraction(sp.together(applyL(CH, S_star)))
chs = [sp.Poly(sp.expand(sp.cancel(gden) * c.as_expr()), t) for c in CH]
ghs = sp.Poly(sp.expand(gnum), t)
mmc, Nmc, _ = bounds(chs, ghs)
solc, msgc = solve_box(chs, ghs, mmc, Nmc)
check("control: manufactured rational S* recovered exactly",
      solc is not None and sp.simplify(solc - S_star) == 0, msgc)
mm0, Nm0, _ = bounds(CH, sp.Poly(0, t))
sol0, msg0 = solve_box(CH, sp.Poly(0, t), mm0, Nm0)
check("control: homogeneous equation has only S = 0 (matches the SL(2,C) lemma)",
      sol0 is not None and sp.simplify(sol0) == 0, msg0)
solf, msgf = solve_box(chs, ghs, {p: 0 for p in mmc}, 3)
check("FALSIFICATION: undersized box does NOT recover S*",
      solf is None or sp.simplify(solf - S_star) != 0, msgf)

# the real run
mm, Nmax, mu = bounds(CH, GH)
for p in sorted(mm, key=lambda z: sp.N(z)):
    nu_p, ind = local_ind(CH, p)
    w = ordp(GH, p)
    emit("    p = %-5s nu_p = %-2d ord_p(g) = %-2s ind. integer roots (-m>0): %-8s m_max = %d"
         % (p, nu_p, (w[0] if w else "-"), pos_roots(ind), mm[p]))
emit("    infinity: mu = %d, numerator degree bound = %d" % (mu, Nmax))
sol, msg = solve_box(CH, GH, mm, Nmax)
if sol is not None:
    residT = sp.simplify(applyL(CH, sol) - GH.as_expr())
    theorem = residT != 0
    emit("    candidate residual: %s" % residT)
else:
    theorem = True
check("THEOREM: Ntilde_v(S2) = 2 Delta has no rational solution", theorem, msg)
emit("")

# =====================================================================
# Stage 5: independence sanity scan (the lemma guarantees all degrees)
# =====================================================================
emit("Stage 5: sanity -- no polynomial relation a W + b W' + c W'' = 0, deg <= 10")
Wp_num = theta(W)[1:] + [F(0)]
W_sh = W[1:] + [F(0)]
Wp = sdiv(Wp_num, W_sh)                       # W' = theta(W)/W, both / q
Wpp = sdiv(theta(Wp)[1:] + [F(0)], W_sh)
cols = []
tp = [F(1)] + [F(0)] * (NQ - 1)
for k in range(11):
    if k:
        tp = smul(tp, tser)
    for base in (W, Wp, Wpp):
        cols.append(smul(tp, base))
nr = min(NQ - 5, len(cols) + 40)
mat = [[cols[c][i] for c in range(len(cols))] for i in range(nr)]
rk = 0
for c in range(len(cols)):
    pz = next((i for i in range(rk, nr) if mat[i][c] != 0), None)
    if pz is None:
        continue
    mat[rk], mat[pz] = mat[pz], mat[rk]
    invv = mat[rk][c]
    mat[rk] = [z / invv for z in mat[rk]]
    for i in range(nr):
        if i != rk and mat[i][c] != 0:
            fz = mat[i][c]
            mat[i] = [mat[i][k] - fz * mat[rk][k] for k in range(len(cols))]
    rk += 1
check("full column rank %d/%d (independence to this degree; all degrees by the lemma)"
      % (rk, len(cols)), rk == len(cols))
emit("")

verdict = ("ALL CHECKS PASS -- Phi = 2 rho1 W/v + 2 Int Delta W/v dt + 2/15 with "
           "Delta = (13t-1)/(30t(t-1)^2), and Phi is NOT in Q(t) + v^{-1} "
           "span{W, W', W''}: the Green's function combination d/dE[(E-1)G_disp] "
           "is weight-two quasimodular, (E-1)G_disp itself is not."
           if ok_all else "CERTIFICATION FAILED -- see FAIL lines above")
emit("RESULT: " + verdict)

with open(os.path.join(HERE, "CERTIFICATE_phi_obstruction.txt"), "w",
          encoding="utf-8") as fh:
    fh.write("\n".join(LINES) + "\n")

sys.exit(0 if ok_all else 1)
