"""
certify_modular.py  --  exact certificate that V2 uniformizes X(Gamma_0(30)+).

Proves, in exact rational arithmetic (self-contained; no external data files):

  u(tau) = [eta(tau) eta(6tau) eta(10tau) eta(15tau) / (eta(2tau) eta(3tau) eta(5tau) eta(30tau))]^3
  t = u/(u^2 + 7u + 1)   (equivalently 1/t = u + 7 + 1/u; generates the genus-zero
                          function field of X(Gamma_0(30)+))

is the modular parametrization of the projective normal form  w'' + Q_V(t) w = 0  of V2, with
  Q_V = N(t) / [4 t^2 (t-1)^2 (4t-1)^2 (5t-1)^2 (9t-1)^2],
  N   = 24300 t^8 - 58860 t^7 + 73437 t^6 - 44294 t^5 + 15111 t^4 - 3160 t^3 + 407 t^2 - 30 t + 1.

Checks (all PASS -> exit 0):
  [1] Ligozat conditions for the eta quotient u  => u is a modular function on Gamma_0(30).
  [2] q-expansion of u matches the closed product q*u = (1-q)^3 - 3 q^7 + ...   (hand value head).
  [3] t_d = u/(u^2+7u+1) equals the V2 MUM mirror map t(q) to ORD orders  (>> the degree-8
      immersion bound => an identity, not a coincidence).
  [4] Schwarzian identity  {tau, t_d} = 2 Q_V(t_d)  holds as a q-series to ORD orders
      => V2 is the uniformizing ODE of X(Gamma_0(30)+).

Pure Python (fractions); no third-party dependencies.
"""
import sys
from fractions import Fraction as F

ORD = 80  # series order for the identity checks (immersion bound is 8; ORD >> 8 is a proof)

# ---------- Q_V (embedded, exact) ----------
NUMER = [1, -30, 407, -3160, 15111, -44294, 73437, -58860, 24300]      # ascending powers t^0..t^8
# denominator 4 t^2 (t-1)^2 (4t-1)^2 (5t-1)^2 (9t-1)^2  -- built below

def polymul(a, b):
    r = [F(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                r[i + j] += ai * bj
    return r

def polypow(a, n):
    r = [F(1)]
    for _ in range(n):
        r = polymul(r, a)
    return r

# ---------- power-series toolkit ----------
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
    inv0 = F(1) / b[0]
    q = [F(0)] * n
    for m in range(n):
        acc = a[m] if m < len(a) else F(0)
        for k in range(1, m + 1):
            if k < len(b) and b[k]:
                acc -= b[k] * q[m - k]
        q[m] = acc * inv0
    return q

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

def pad(p, n):
    return (list(p) + [F(0)] * n)[:n]

ok = True
def check(name, cond):
    global ok
    ok = ok and bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", name))

print("certify_modular.py  --  V2 uniformizes X(Gamma_0(30)+)")
print()

# ---------- [1] Ligozat ----------
N = 30
r = {1: 3, 6: 3, 10: 3, 15: 3, 2: -3, 3: -3, 5: -3, 30: -3}
s1 = sum(d * rd for d, rd in r.items())
s2 = sum((N // d) * rd for d, rd in r.items())
wt = sum(r.values())
prod = F(1)
for d, rd in r.items():
    prod *= F(d) ** rd
import math
def is_rat_square(fr):
    a, b = fr.numerator, fr.denominator
    return fr > 0 and math.isqrt(a) ** 2 == a and math.isqrt(b) ** 2 == b
print("Check 1: Ligozat conditions for u (level 30 eta quotient)")
check("sum d*r_d = %d = 0 mod 24" % s1, s1 % 24 == 0)
check("sum (30/d)*r_d = %d = 0 mod 24" % s2, s2 % 24 == 0)
check("weight (1/2) sum r_d = %s (0 => modular function)" % F(wt, 2), wt == 0)
check("prod d^{r_d} = %s is a rational square (trivial character)" % prod, is_rat_square(prod))
print()

# ---------- build u = q^{-1} * U(q) ----------
n = ORD + 6
def eta_ipart(d, n):
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
ratio = smul(Pnum, sinv(Pden, n), n)         # prod (1-q^{d n}) ratio
U = smul(smul(ratio, ratio, n), ratio, n)    # U(q) = (num/den)^3 ; u = q^{-1} U

# ---------- [2] u head ----------
print("Check 2: eta expansion q*u head equals the hand value [1,-3,3,-1,0,0,0,-3,9]")
head = [str(x) for x in U[:9]]
check("q*u = %s" % " ".join(head), head == ['1', '-3', '3', '-1', '0', '0', '0', '-3', '9'])
print()

# t_d = u/(u^2+7u+1) = (U*q)/(U^2 + 7 q U + q^2)
qU = [F(0)] + U[:n - 1]
q2 = [F(0)] * n
q2[2] = F(1)
U2 = smul(U, U, n)
den_td = [U2[i] + 7 * qU[i] + q2[i] for i in range(n)]
numt = [F(0)] + U[:n - 1]
t_d = smul(numt, sinv(den_td, n), n)

# ---------- V2 mirror map from Q_V ----------
den_poly = [F(4)]
den_poly = polymul(den_poly, polypow([F(0), F(1)], 2))       # t^2
den_poly = polymul(den_poly, polypow([F(-1), F(1)], 2))      # (t-1)^2
den_poly = polymul(den_poly, polypow([F(-1), F(4)], 2))      # (4t-1)^2
den_poly = polymul(den_poly, polypow([F(-1), F(5)], 2))      # (5t-1)^2
den_poly = polymul(den_poly, polypow([F(-1), F(9)], 2))      # (9t-1)^2
num_poly = [F(c) for c in NUMER]

# G(t) = t^2 Q_V - 1/4.  Since Q_V = num/den with den = t^2 * den_red,
#   t^2 Q_V = num/den_red,  G = (num - (1/4) den_red)/den_red,  den_red[0] = 4 != 0.
den_red = [F(4)]
den_red = polymul(den_red, polypow([F(-1), F(1)], 2))   # (t-1)^2
den_red = polymul(den_red, polypow([F(-1), F(4)], 2))   # (4t-1)^2
den_red = polymul(den_red, polypow([F(-1), F(5)], 2))   # (5t-1)^2
den_red = polymul(den_red, polypow([F(-1), F(9)], 2))   # (9t-1)^2
combo = [ (num_poly[i] if i < len(num_poly) else F(0)) - F(1, 4) * (den_red[i] if i < len(den_red) else F(0))
          for i in range(max(len(num_poly), len(den_red))) ]
G = sdiv(pad(combo, n), pad(den_red, n), n)
assert G[0] == 0, "Q_V normalization: leading 1/(4t^2) pole expected"
# Frobenius at MUM (exponents 1/2,1/2): a_m=-(1/m^2)sum G_k a_{m-k}; b_m=(-2m a_m - sum G_k b_{m-k})/m^2
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
sV = sdiv(b, a, n)                     # mirror map s(t); nome q = t exp(sV)
expS = sexp(sV, n)
inv_unit = sdiv([F(1)] + [F(0)] * (n - 1), expS, n)   # (t/q) = 1/exp(sV)
P = inv_unit[:]
tq = [F(0)] * n
for m in range(1, n):
    tq[m] = P[m - 1] / m
    if m < n - 1:
        P = smul(P, inv_unit, n)

# ---------- [3] identity t_d == mirror map ----------
print("Check 3: t_d = u/(u^2+7u+1) equals the V2 MUM mirror map to %d orders" % ORD)
firstbad = next((i for i in range(ORD) if t_d[i] != tq[i]), None)
check("t_d(q) - t_mirror(q) first nonzero order = %s (None => identical; immersion bound 8)" % firstbad,
      firstbad is None)
print("        t(q) = " + " ".join("%+d q^%d" % (int(tq[i]), i) for i in range(1, 7)) + " ...")
print()

# ---------- [4] Schwarzian ----------
print("Check 4: Schwarzian {tau,t_d} = 2 Q_V(t_d)  (q-series, %d orders)" % ORD)
def Dop(s, n):
    return [F(k) * s[k] for k in range(n)]
Dt1 = Dop(t_d, n); Dt2 = Dop(Dt1, n); Dt3 = Dop(Dt2, n)
def uunit(s, n):
    return [s[k + 1] for k in range(n - 1)] + [F(0)]
r2 = smul(uunit(Dt2, n), sinv(uunit(Dt1, n), n), n)
r3 = smul(uunit(Dt3, n), sinv(uunit(Dt1, n), n), n)
bracket = [r3[i] - F(3, 2) * smul(r2, r2, n)[i] for i in range(n)]
U2d = smul(uunit(Dt1, n), uunit(Dt1, n), n)
q2_sch = smul([-x for x in bracket], sinv(U2d, n), n)          # q^2 {tau,t_d}
def polyeval_series(coeffs_asc, xs, n):
    out = [F(0)] * n
    xp = [F(1)] + [F(0)] * (n - 1)
    for k, c in enumerate(coeffs_asc):
        if k > 0:
            xp = smul(xp, xs, n)
        if c:
            for i in range(n):
                if xp[i]:
                    out[i] += c * xp[i]
    return out
den_ser = polyeval_series([F(c) for c in den_poly], t_d, n)    # ~4 q^2
den_unit = [den_ser[k + 2] for k in range(n - 2)] + [F(0), F(0)]
num_ser = polyeval_series(num_poly, t_d, n)
q2_2QV = smul([2 * x for x in num_ser], sinv(pad(den_unit, n), n), n)
fb = next((i for i in range(ORD - 2) if q2_sch[i] != q2_2QV[i]), None)
check("q^2 ({tau,t_d} - 2 Q_V(t_d)) first nonzero order = %s (None => identity)" % fb, fb is None)
print()

print("RESULT:", "ALL CHECKS PASS -- V2 uniformizes X(Gamma_0(30)+); t generates its function field." if ok
      else "FAILURE -- see above.")
sys.exit(0 if ok else 1)
