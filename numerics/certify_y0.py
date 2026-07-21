"""
certify_y0.py -- exact certificate for the closed form of the weight-two period
y0 = Phi'/2 of the hyperkagome lattice Green's function (level-30 modular form).

THEOREM (proven; this script packages the exact checks).  With q the nome and
  u   = [eta(tau)eta(6tau)eta(10tau)eta(15tau)/(eta(2tau)eta(3tau)eta(5tau)eta(30tau))]^3,
  ell = q dlog u
      = (1/8)[E2(q)+6E2(q^6)+10E2(q^10)+15E2(q^15)-2E2(q^2)-3E2(q^3)-5E2(q^5)-30E2(q^30)],
  t   = u/(u^2+7u+1)      (Hauptmodul of Gamma_0(30)+; equals the MUM mirror map of Q_V),
  W   = q dt/dq,
  v   = sqrt((1-4t)(1-5t)(1-9t)),   v(0) = +1,
the weight-two period  y0(t) = 1 + 10t + 87t^2 + 724t^3 + ...  (the holomorphic solution
at t=0 of the order-3 Picard-Fuchs operator M of the lattice Green's function) satisfies

      B(t) * (y0 o t)(q) * v * W  =  A0(t) * W^2 + A1(t) * (q dW/dq),

  B  = -30 t^2 (t-1)^2 (4t-1)(5t-1)(9t-1),
  A0 = 4050 t^6 + 2445 t^5 - 11436 t^4 + 8000 t^3 - 2130 t^2 + 231 t - 8,
  A1 = -t (t-1)(4t-1)(5t-1)(9t-1)(15 t^2 + 17 t - 8);

equivalently  y0 = [rho0 W + rho1 W'] / v,  W' = dW/dt = (q dW/dq)/W,
  rho1 = (15 t^2 + 17 t - 8) / (30 t (t-1)),
  rho0 = -A0 / (30 t^2 (t-1)^2 (4t-1)(5t-1)(9t-1)).

PROOF CHAIN:
  [I]   Schwarzian identity {tau,t} = 2 Q_V (certify_modular.py) => w1 = (dt/dtau)^{1/2}
        solves the projective normal form  w'' + Q_V w = 0.
  [II]  Sym^2 lemma: y'' = -Q y => (y^2)''' + 4 Q (y^2)' + 2 Q' (y^2) = 0 => Ntilde(W) = 0,
        Ntilde = D^3 + 4 Q_V D + 2 Q_V'.
  [III] M_v := v o M o v^{-1} = sum_i c_i (D - w)^i,  w = (1/2) dlog((1-4t)(1-5t)(1-9t)),
        has coefficients in Q(t) and annihilates y0*v.
  [IV]  Ore identity  M_v . T = X . Ntilde  in Q(t)<D>, remainder of the right division
        exactly 0, T = rho0 + rho1 D  =>  T maps Sol(Ntilde) -> Sol(M_v), so T(W) solves M_v.
  [V]   MUM at t=0 (CERTIFICATE_mum.txt): the analytic solution line is 1-dimensional;
        y0*v and T(W) are both analytic at q=0 with equal q^0 coefficient => y0*v = T(W). QED.

Checks (all PASS -> exit 0):
  [1] y0 (nu.json head, extended by M's own recurrence) satisfies M(y0) = 0 exactly.
  [2] eta<->E2: the ell combination equals q dlog u to ORD orders.
  [3] t = u/(u^2+7u+1) equals the MUM mirror map of Q_V to ORD orders.
  [4] v^2 = (1-4t)(1-5t)(1-9t) composed with t(q); v(0) = +1.
  [5] MAIN IDENTITY  B (y0 o t) v W = A0 W^2 + A1 (q dW/dq)  to ORD orders.
  [6] (sympy) the Ore operator identity  M_v . T - X . Ntilde = 0  exactly.

Self-contained: pure-stdlib Fraction power series (never floats); loads only nu.json and
M_coeffs.json from its own directory; sympy is used ONLY for check [6].  ASCII output only.
"""
import sys, os, json
from fractions import Fraction as F

ORD = 120          # series depth for the q-series identity checks (immersion/fit bounds << ORD)
NY  = 135          # number of t-series terms of y0 (nu.json head + M-recurrence extension)
NQ  = 130          # number of q-series terms

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------- polynomial toolkit (ascending Fraction coeff lists) ----------------
def polymul(a, b):
    r = [F(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    r[i + j] += ai * bj
    return r

def polypow(a, n):
    r = [F(1)]
    for _ in range(n):
        r = polymul(r, a)
    return r

def pad(p, n):
    return (list(p) + [F(0)] * n)[:n]

# ---------------- power-series toolkit (exact, ascending in q) ----------------
def smul(a, b, n):
    r = [F(0)] * n
    for i, ai in enumerate(a[:n]):
        if ai:
            for j, bj in enumerate(b[:n - i]):
                if bj:
                    r[i + j] += ai * bj
    return r

def sinv(a, n):
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

def sdiv(a, b, n):
    return smul(a, sinv(b, n), n)

def sexp(s, n):
    u = [F(0)] * n
    u[0] = F(1)
    for m in range(1, n):
        acc = F(0)
        for k in range(1, m + 1):
            if s[k]:
                acc += k * s[k] * u[m - k]
        u[m] = acc / m
    return u

def ssqrt(a, n):
    """principal square root of a series with a[0] = 1."""
    assert a[0] == 1
    s = [F(0)] * n
    s[0] = F(1)
    for m in range(1, n):
        acc = a[m] if m < len(a) else F(0)
        for k in range(1, m):
            acc -= s[k] * s[m - k]
        s[m] = acc / 2
    return s

def theta(s, n):
    """q d/dq acting on a series (theta operator)."""
    return [F(k) * s[k] for k in range(n)]

def scomp(f, g, n):
    """f(g(q)) for polynomial/series f (in t) and series g with g[0] = 0 (Horner)."""
    assert g[0] == 0
    acc = [F(0)] * n
    for k in range(min(len(f), n) - 1, -1, -1):
        acc = smul(acc, g, n)
        acc[0] += f[k]
    return acc

def polyser(ci, ts, n):
    """evaluate a polynomial (Fraction coeff list, ascending) at the series ts(q)."""
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

# ---------------- PASS/FAIL harness ----------------
ok = True
def check(name, cond):
    global ok
    ok = ok and bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", name))

print("certify_y0.py  --  closed form of the weight-two period y0 = Phi'/2")
print()

# ================= operator M (level-30 Picard-Fuchs), from M_coeffs.json =================
Mj = json.load(open(os.path.join(HERE, "M_coeffs.json")))
assert Mj["order"] == 3, "M_coeffs.json: expected order 3"
Mc = [[F(int(x)) for x in lst] for lst in Mj["coeffs"]]     # c[0..3] ascending in t

def series_deriv(a, k=1):
    r = list(a)
    for _ in range(k):
        r = [F(m + 1) * r[m + 1] for m in range(len(r) - 1)] if len(r) > 1 else [F(0)]
    return r

def apply_M_series(y, c=Mc):
    """coefficient list of M(y) for a power series y (M = sum_i c_i(t) d^i/dt^i)."""
    N = len(y)
    acc = [F(0)] * (N + 20)
    for i in range(4):
        di = series_deriv(y, i)
        for m, cm in enumerate(c[i]):
            if cm:
                for k, dk in enumerate(di):
                    if m + k < len(acc):
                        acc[m + k] += cm * dk
    return acc

def y0_extend(y0_known, N, c=Mc):
    """extend the holomorphic solution y0 to N terms via M(y0) = 0 as a recurrence.
    At t=0 the exponents are {-1,-1,0}; y0 is the exponent-0 (holomorphic) solution, so its
    coefficient a[n0] is fixed by the equation 'coeff of t^{n0-1} in M(y) = 0'.  Because
    c3 starts at t^2, c2 at t^1, c1 at t^0, no a_n with n > n0 enters that row -- march upward."""
    a = list(y0_known) + [None] * (N - len(y0_known))
    def ff(n, i):
        r = F(1)
        for k in range(i):
            r *= (n - k)
        return r
    ci = [dict(enumerate(ci_)) for ci_ in c]
    for n0 in range(len(y0_known), N):
        m = n0 - 1
        coef_unknown = F(0)
        rhs = F(0)
        for i in range(4):
            for j, cij in ci[i].items():
                if cij == 0:
                    continue
                n = m + i - j
                if n < 0:
                    continue
                val = cij * ff(n, i)
                if n == n0:
                    coef_unknown += val
                elif n < n0:
                    if a[n] is None:
                        raise RuntimeError("need a[%d] at n0=%d" % (n, n0))
                    rhs += val * a[n]
                else:
                    raise RuntimeError("future coeff a[%d] at n0=%d" % (n, n0))
        if coef_unknown == 0:
            raise RuntimeError("singular recurrence at n0=%d" % n0)
        a[n0] = -rhs / coef_unknown
    return a

# ---------------- Check 1 : y0 head + M-recurrence, residual exactly 0 ----------------
print("Check 1: y0 (nu.json head + M-recurrence) satisfies M(y0) = 0 exactly")
nu = [F(p, qd) for p, qd in json.load(open(os.path.join(HERE, "nu.json")))]
y0_head = [F(k + 1) * nu[k + 1] / 2 for k in range(len(nu) - 1)]      # y0[k] = (k+1) nu[k+1]/2
check("y0 head = (k+1) nu[k+1]/2 = [1,10,87,724]", [str(x) for x in y0_head[:4]] == ['1', '10', '87', '724'])
y0 = y0_extend(y0_head, NY)
check("y0 extended to %d terms; every coefficient integral" % NY, all(x.denominator == 1 for x in y0))
res = apply_M_series(y0)
bad1 = next((m for m in range(ORD + 6) if res[m] != 0), None)
check("M(y0) = 0 through t^%d (first nonzero order = %s)" % (ORD + 5, bad1), bad1 is None)
print()

# ================= q-series ingredients: u, ell, t, W, v =================
n = NQ
def eta_ipart(d, n):
    """integer part of prod_{k>=1}(1 - q^{d k}) as a q-series to n terms."""
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

Pnum = [F(1)] + [F(0)] * (n - 1)
for d in (1, 6, 10, 15):
    Pnum = smul(Pnum, eta_ipart(d, n), n)
Pden = [F(1)] + [F(0)] * (n - 1)
for d in (2, 3, 5, 30):
    Pden = smul(Pden, eta_ipart(d, n), n)
ratio = sdiv(Pnum, Pden, n)
U = smul(smul(ratio, ratio, n), ratio, n)                  # u = U/q  (so q*u = U)

# ---------------- Check 2 : eta <-> E2  (ell = q dlog u) ----------------
print("Check 2: eta<->E2 identity  ell = q dlog u  to %d orders" % ORD)
check("q*u head = [1,-3,3,-1,0,0,0,-3,9]",
      [str(x) for x in U[:9]] == ['1', '-3', '3', '-1', '0', '0', '0', '-3', '9'])
def sigma1(m):
    return sum(dd for dd in range(1, m + 1) if m % dd == 0)
def E2(dd, n):
    e = [F(0)] * n
    e[0] = F(1)
    for m in range(1, n):
        if m % dd == 0:
            e[m] = F(-24) * sigma1(m // dd)
    return e
ell = [F(0)] * n
for dd, sgn in ((1, 1), (6, 1), (10, 1), (15, 1), (2, -1), (3, -1), (5, -1), (30, -1)):
    Ed = E2(dd, n)
    for i in range(n):
        ell[i] += sgn * F(dd, 8) * Ed[i]
qdlU = sdiv(theta(U, n), U, n)                             # q dlog U
ell_direct = [qdlU[i] - (F(1) if i == 0 else F(0)) for i in range(n)]   # q dlog u = q dlog U - 1
bad2 = next((i for i in range(ORD) if ell[i] != ell_direct[i]), None)
check("(1/8)[E2(q)+6E2(q^6)+10E2(q^10)+15E2(q^15)-2E2(q^2)-3E2(q^3)-5E2(q^5)-30E2(q^30)]"
      " = q dlog u (first bad order = %s)" % bad2, bad2 is None)
print()

# t = u/(u^2+7u+1) = qU/(U^2 + 7 q U + q^2)   [eta route]
qU = [F(0)] + U[:n - 1]
q2 = [F(0)] * n
q2[2] = F(1)
U2 = smul(U, U, n)
den_td = [U2[i] + 7 * qU[i] + q2[i] for i in range(n)]
t_eta = sdiv(qU, den_td, n)

# ---------------- Check 3 : t equals the MUM mirror map of Q_V ----------------
print("Check 3: t = u/(u^2+7u+1) equals the MUM mirror map of Q_V to %d orders" % ORD)
NUMER = [1, -30, 407, -3160, 15111, -44294, 73437, -58860, 24300]     # Q_V numerator (embedded)
den_red = [F(4)]
for (c0, c1) in ((-1, 1), (-1, 4), (-1, 5), (-1, 9)):                 # 4 (t-1)^2 (4t-1)^2 (5t-1)^2 (9t-1)^2
    den_red = polymul(den_red, polymul([F(c0), F(c1)], [F(c0), F(c1)]))
combo = [(F(NUMER[i]) if i < len(NUMER) else F(0)) - F(1, 4) * (den_red[i] if i < len(den_red) else F(0))
         for i in range(max(len(NUMER), len(den_red)))]
G = sdiv(pad(combo, n), pad(den_red, n), n)                # G = t^2 Q_V - 1/4  (regular at 0)
assert G[0] == 0, "Q_V normalization: expected leading 1/(4 t^2) pole"
# Frobenius at MUM (double exponent 1/2): a_m = -(1/m^2) sum G_k a_{m-k};
#                                          b_m = (-2 m a_m - sum G_k b_{m-k})/m^2
a = [F(0)] * n
b = [F(0)] * n
a[0] = F(1)
for m in range(1, n):
    aa = F(0)
    bb = F(0)
    for k in range(1, m + 1):
        if G[k]:
            aa += G[k] * a[m - k]
            bb += G[k] * b[m - k]
    a[m] = -aa / (m * m)
    b[m] = (F(-2 * m) * a[m] - bb) / (m * m)
sV = sdiv(b, a, n)                                         # mirror map s(t); nome q = t exp(sV)
expS = sexp(sV, n)
inv_unit = sdiv([F(1)] + [F(0)] * (n - 1), expS, n)        # (t/q) = 1/exp(sV)
P = inv_unit[:]
tq = [F(0)] * n
for m in range(1, n):
    tq[m] = P[m - 1] / m
    if m < n - 1:
        P = smul(P, inv_unit, n)
bad3 = next((i for i in range(ORD) if t_eta[i] != tq[i]), None)
check("t_eta(q) - t_mirror(q) first nonzero order = %s (None => identity; immersion bound 8)" % bad3,
      bad3 is None)
print("        t(q) = " + " ".join("%+d q^%d" % (int(t_eta[i]), i) for i in range(1, 7)) + " ...")
print()

t_ser = t_eta      # canonical t(q) for the downstream identity

# W = q dt/dq ;  v = sqrt((1-4t)(1-5t)(1-9t)) o t(q)
W = theta(t_ser, n)
prod_poly = polymul(polymul([F(1), F(-4)], [F(1), F(-5)]), [F(1), F(-9)])   # (1-4t)(1-5t)(1-9t)
prod_q = scomp(prod_poly, t_ser, n)
V = ssqrt(prod_q, n)

# ---------------- Check 4 : v^2 identity, v(0) = +1 ----------------
print("Check 4: v^2 = (1-4t)(1-5t)(1-9t) composed with t(q); v(0) = +1")
vv = smul(V, V, n)
check("v(0) = +1", V[0] == 1)
bad4 = next((i for i in range(ORD) if vv[i] != prod_q[i]), None)
check("v^2 - (1-4t)(1-5t)(1-9t) o t(q) first nonzero order = %s" % bad4, bad4 is None)
print()

# ================= closed-form data B, A0, A1 (factored theorem forms) =================
B_INT  = [0, 0, 30, -600, 4140, -12000, 13830, -5400]                 # expanded literals (guards)
A0_INT = [-8, 231, -2130, 8000, -11436, 2445, 4050]
A1_INT = [0, 8, -169, 1260, -3986, 4432, 1155, -2700]
# B = -30 t^2 (t-1)^2 (4t-1)(5t-1)(9t-1)
Bpoly = polymul(polypow([F(0), F(1)], 2), polypow([F(-1), F(1)], 2))
for fac in ([F(-1), F(4)], [F(-1), F(5)], [F(-1), F(9)]):
    Bpoly = polymul(Bpoly, fac)
Bpoly = [F(-30) * x for x in Bpoly]
assert [int(x) for x in Bpoly] == B_INT, "B factored form != expanded literal"
# A1 = -t (t-1)(4t-1)(5t-1)(9t-1)(15 t^2 + 17 t - 8)
A1poly = [F(0), F(1)]
for fac in ([F(-1), F(1)], [F(-1), F(4)], [F(-1), F(5)], [F(-1), F(9)], [F(-8), F(17), F(15)]):
    A1poly = polymul(A1poly, fac)
A1poly = [F(-1) * x for x in A1poly]
assert [int(x) for x in A1poly] == A1_INT, "A1 factored form != expanded literal"
A0poly = [F(c) for c in A0_INT]                            # A0 is irreducible over Q -- explicit sextic

# ---------------- Check 5 : main identity ----------------
print("Check 5: MAIN IDENTITY  B*(y0 o t)*v*W = A0*W^2 + A1*(q dW/dq)  to %d orders" % ORD)
Y = scomp(y0, t_ser, n)                                    # y0(t(q))
check("y0(t(q)) head = [1,10,47,148,407,1254]",
      [str(x) for x in Y[:6]] == ['1', '10', '47', '148', '407', '1254'])
B_ser  = polyser(Bpoly,  t_ser, n)
A0_ser = polyser(A0poly, t_ser, n)
A1_ser = polyser(A1poly, t_ser, n)
lhs = smul(smul(B_ser, smul(Y, V, n), n), W, n)
thW = theta(W, n)                                          # q dW/dq
rhs = [smul(A0_ser, smul(W, W, n), n)[i] + smul(A1_ser, thW, n)[i] for i in range(n)]
bad5 = next((i for i in range(ORD) if lhs[i] != rhs[i]), None)
check("B*(y0 o t)*v*W - (A0*W^2 + A1*(q dW/dq)) first nonzero order = %s" % bad5, bad5 is None)
FREE = 32     # B,A0,A1,A2 ansatz: 4 polynomials of degree 7 -> 4*8 = 32 free coefficients
print("        margin = %d orders checked - %d free coefficients = %d  (overdetermined => identity)"
      % (ORD, FREE, ORD - FREE))
print()

# ---------------- Check 6 : symbolic Ore operator identity (sympy) ----------------
print("Check 6: Ore operator identity  M_v . T = X . Ntilde  (exact rational functions)")
import sympy as sp
t = sp.symbols('t')

def opmul(A, B):
    """product in Q(t)<D>:  op = [a0,a1,...] = sum ai D^i.  D^i o (bj D^j) = sum C(i,k) bj^(k) D^{i-k+j}."""
    res = [sp.Integer(0)] * (len(A) + len(B) - 1)
    for i, ai in enumerate(A):
        if ai == 0:
            continue
        for j, bj in enumerate(B):
            if bj == 0:
                continue
            for k in range(i + 1):
                term = sp.binomial(i, k) * sp.diff(bj, t, k)
                if term != 0:
                    res[i - k + j] = sp.cancel(res[i - k + j] + ai * term)
    return [sp.cancel(x) for x in res]

def opadd(A, B):
    m = max(len(A), len(B))
    return [sp.cancel((A[i] if i < len(A) else 0) + (B[i] if i < len(B) else 0)) for i in range(m)]

def opscale(A, s):
    return [sp.cancel(s * x) for x in A]

def optrim(A):
    A = list(A)
    while len(A) > 1 and A[-1] == 0:
        A.pop()
    return A

def rdiv(Pop, Nop):
    """right-divide Pop by leading-unit Nop:  Pop = X.Nop + R,  deg R < deg Nop."""
    Pop = [sp.cancel(x) for x in Pop]
    dN = len(Nop) - 1
    lcN = Nop[-1]
    X = [sp.Integer(0)] * max(1, len(Pop) - dN)
    while len(optrim(Pop)) - 1 >= dN:
        Pop = optrim(Pop)
        m = len(Pop) - 1
        coef = sp.cancel(Pop[m] / lcN)
        X[m - dN] = sp.cancel(X[m - dN] + coef)
        mono = [sp.Integer(0)] * (m - dN) + [coef]
        Pop = optrim(opadd(Pop, opscale(opmul(mono, Nop), -1)))
        if all(x == 0 for x in Pop):
            break
    return X, Pop

# M (from M_coeffs.json);  M_v = v o M o v^{-1} = sum_i c_i (D - w)^i,  w = v'/v
c = [sum(sp.Integer(int(x)) * t**k for k, x in enumerate(lst)) for lst in Mj["coeffs"]]
w = sp.cancel(sp.Rational(1, 2) * sp.diff(sp.log((1 - 4*t) * (1 - 5*t) * (1 - 9*t)), t))
Dw = [-w, sp.Integer(1)]                                   # D - w
pw = [[sp.Integer(1)]]
for i in range(3):
    pw.append(opmul(Dw, pw[-1]))
Mv = [sp.Integer(0)]
for i in range(4):
    Mv = opadd(Mv, opscale(pw[i], c[i]))
Mv = optrim(Mv)
assert len(Mv) - 1 == 3, "M_v is not order 3"

# T = rho0 + rho1 D ;  rho0 = A0/B,  rho1 = A1/B ;  Ntilde = [2 Q_V', 4 Q_V, 0, 1]
Bpol_s = sum(sp.Integer(cc) * t**k for k, cc in enumerate(B_INT))
A0_s   = sum(sp.Integer(cc) * t**k for k, cc in enumerate(A0_INT))
A1_s   = sum(sp.Integer(cc) * t**k for k, cc in enumerate(A1_INT))
rho0 = sp.cancel(A0_s / Bpol_s)
rho1 = sp.cancel(A1_s / Bpol_s)
T = [rho0, rho1]
QV = sp.cancel(sum(sp.Integer(cc) * t**k for k, cc in enumerate(NUMER)) /
               (4 * t**2 * (t - 1)**2 * (4*t - 1)**2 * (5*t - 1)**2 * (9*t - 1)**2))
Ntil = [sp.cancel(2 * sp.diff(QV, t)), sp.cancel(4 * QV), sp.Integer(0), sp.Integer(1)]

MvT = opmul(Mv, T)
X, R = rdiv(MvT, Ntil)
R = optrim(R)
allzero = all(sp.simplify(x) == 0 for x in R)
print("  rho1 =", rho1)
print("  rho0 =", rho0)
print("  order(M_v) =", len(optrim(Mv)) - 1, "; order(M_v . T) =", len(optrim(MvT)) - 1,
      "; order(Ntilde) =", len(optrim(Ntil)) - 1)
check("right-division remainder R = (M_v . T) mod Ntilde is exactly 0 (entries %s)"
      % [sp.simplify(x) for x in R], allzero)
for i, xi in enumerate(optrim(X)):
    print("    X[%d] = %s" % (i, sp.factor(sp.cancel(xi))))
print("  [V] MUM anchor: (y0*v)(q=0) =", Y[0] * V[0],
      "; analytic line 1-dim (CERTIFICATE_mum.txt) => y0*v = T(W)")
print()

print("RESULT:", "ALL CHECKS PASS -- y0 = [rho0 W + rho1 W']/v ;  B y0 v W = A0 W^2 + A1 (q dW/dq)." if ok
      else "FAILURE -- see above.")
sys.exit(0 if ok else 1)
