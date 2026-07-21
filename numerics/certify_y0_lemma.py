#!/usr/bin/env python3
"""Certificate for the exclusion lemma behind remark (ii) of the y0 section:

    y0 is not (meromorphic modular form of any weight, any finite-order
    multiplier, any level commensurable with Gamma_0(30)+) x (algebraic
    function); in particular y0 is not a pure eta quotient.

The proof in the paper reduces to two computational facts certified here.

CHECK 1 (load-bearing).  The order-one intertwiner T = rho0 + rho1*d/dt of the
  proven Ore identity M_v T = X Ntilde is INJECTIVE on Sol(Ntilde), where
  Ntilde = D^3 + 4*Q_V*D + 2*Q_V' is the symmetric square of w'' + Q_V w = 0.
  Any kernel element is a multiple of z = exp(-int rho0/rho1), and z solves
  Ntilde iff E(t) := Ntilde(z)/z = r'' + 3rr' + r^3 + 4*Q_V*r + 2*Q_V'
  (r = -rho0/rho1 = z'/z) vanishes identically.  We compute, exactly,

      E(t) = -30 * p7(t) / [ t^2 (t-1)(4t-1)(5t-1)(9t-1)(15t^2+17t-8)^3 ],

  with p7 the degree-seven apparent-singularity polynomial of M (the same
  polynomial certified apparent in certify_p7_apparent.py, and the same locus
  carried by the cofactor X of the Ore identity).  E != 0, so T is a bijection
  Sol(Ntilde) -> Sol(M_v): the Galois representation on Sol(M_v) is the
  symmetric square of the standard representation of Gal(w''+Q_V w = 0),
  which together with Gal(M)^0 = SO(3,C) forces Gal(w''+Q_V w=0) = SL(2,C).

CHECK 2 (corroborating scan).  The theorem implies L = W'/W (W = q dt/dq,
  ' = d/dt) is transcendental over C(t).  We certify that no nonzero
  P in Q[t,X] with deg_t <= 28, deg_X <= 12 satisfies P(t, t*L) = 0, by
  exhibiting full column rank (modulo the prime 2^61-1) of the coefficient
  linear system against 480 exact q-expansion orders.  Full column rank mod p
  implies full column rank over Q (any putative integer relation, scaled to
  content one, would survive reduction mod p).

CHECK 0 (pipeline guards).  t(q) built from the level-30 eta quotient
  reproduces the certified Hauptmodul heads, and Ntilde(W) = 0 holds through
  q^553 in the pole-free polynomial form D^2 W''' + 4NDW' + 2(N'D-ND')W = 0
  -- tying the Q_V used here to the modular data end to end.

Run:  python certify_y0_lemma.py      (exit 0 on success; ~2-3 min)
"""
import sys
import sympy as sp

t = sp.Symbol('t')
fails = []
def check(name, ok, detail=""):
    print(("PASS" if ok else "FAIL"), "-", name, ("| " + detail if detail else ""))
    if not ok:
        fails.append(name)

# ---- exact data (as printed in the paper) -----------------------------------
N_QV = (24300*t**8 - 58860*t**7 + 73437*t**6 - 44294*t**5 + 15111*t**4
        - 3160*t**3 + 407*t**2 - 30*t + 1)
D_QV = 4*t**2*(t-1)**2*(4*t-1)**2*(5*t-1)**2*(9*t-1)**2
QV = N_QV / D_QV
rho1 = (15*t**2 + 17*t - 8) / (30*t*(t-1))
rho0 = -(4050*t**6 + 2445*t**5 - 11436*t**4 + 8000*t**3 - 2130*t**2 + 231*t - 8) / \
        (30*t**2*(t-1)**2*(4*t-1)*(5*t-1)*(9*t-1))
# p7 (ascending), identical to certify_p7_apparent.py
P7_ASC = [-64, 1800, -17580, 93657, -290956, 455798, -369600, 101025]
p7 = sum(c*t**i for i, c in enumerate(P7_ASC))

# ---------------------------------------------------- CHECK 1: E(t) explicitly
r = sp.together(-rho0/rho1)
E = sp.together(sp.diff(r, t, 2) + 3*r*sp.diff(r, t) + r**3 + 4*QV*r + 2*sp.diff(QV, t))
E_closed = -30*p7 / (t**2*(t-1)*(4*t-1)*(5*t-1)*(9*t-1)*(15*t**2+17*t-8)**3)
check("E(t) = -30*p7 / [t^2(t-1)(4t-1)(5t-1)(9t-1)(15t^2+17t-8)^3] exactly",
      sp.simplify(sp.cancel(E - E_closed)) == 0)
check("E(t) != 0  (kernel line of T does not solve Ntilde => T injective)",
      sp.cancel(E) != 0, "E(1/2) = %s" % sp.nsimplify(E.subs(t, sp.Rational(1, 2))))

# ---- exact q-series machinery mod a 61-bit prime ----------------------------
P = (1 << 61) - 1
NORD = 560

def euler(nord):
    c = [0]*nord; c[0] = 1
    k = 1
    while True:
        g1 = k*(3*k-1)//2; g2 = k*(3*k+1)//2
        if g1 >= nord and g2 >= nord:
            break
        s = -1 if k % 2 else 1
        if g1 < nord: c[g1] += s
        if g2 < nord: c[g2] += s
        k += 1
    return c

def smul(a, b, nord=NORD):
    out = [0]*nord
    for i, ai in enumerate(a):
        if ai:
            lim = nord - i
            for j in range(min(len(b), lim)):
                if b[j]:
                    out[i+j] = (out[i+j] + ai*b[j]) % P
    return out

def sinv(a, nord=NORD):
    inv0 = pow(a[0], P-2, P)
    out = [0]*nord; out[0] = inv0
    for n in range(1, nord):
        s = 0
        for k in range(1, min(n, len(a)-1) + 1):
            s = (s + a[k]*out[n-k]) % P
        out[n] = (-inv0*s) % P
    return out

def spow(a, e, nord=NORD):
    out = [0]*nord; out[0] = 1
    base = a[:]
    while e:
        if e & 1: out = smul(out, base, nord)
        e >>= 1
        if e: base = smul(base, base, nord)
    return out

def theta(a):
    return [(i*ai) % P for i, ai in enumerate(a)]

def scale(a, d, nord=NORD):
    out = [0]*nord
    for i, ai in enumerate(a):
        if i*d < nord: out[i*d] = ai % P
    return out

phi = euler(NORD)
num = smul(smul(scale(phi, 1), scale(phi, 6)), smul(scale(phi, 10), scale(phi, 15)))
den = smul(smul(scale(phi, 2), scale(phi, 3)), smul(scale(phi, 5), scale(phi, 30)))
qu = spow(smul(num, sinv(den)), 3)                      # q*u, u the eta quotient

qu2 = smul(qu, qu)                                      # t = u/(u^2+7u+1)
den_t = [(qu2[i] + 7*(qu[i-1] if i >= 1 else 0) + (1 if i == 2 else 0)) % P
         for i in range(NORD)]
num_t = [0]*NORD
for i in range(NORD-1): num_t[i+1] = qu[i]
tser = smul(num_t, sinv(den_t))

check("t(q) heads = q - 4q^2 + 12q^3 - 34q^4 + 90q^5 (eta pipeline)",
      [tser[i] % P for i in range(1, 6)] == [x % P for x in (1, -4, 12, -34, 90)])

W = theta(tser)                                          # W = q dt/dq
Wsh = [W[i+1] for i in range(NORD-1)] + [0]              # W/q, a unit
Winv = sinv(Wsh)
def divW(a):
    ash = [a[i+1] for i in range(NORD-1)] + [0]
    return smul(ash, Winv)
def ddt(a):
    return divW(theta(a))

# ------------------------- CHECK 0: Ntilde(W) = 0, pole-free polynomial form
def polyser(expr):
    coeffs = sp.Poly(sp.expand(expr), t).all_coeffs()[::-1]
    out = [0]*NORD
    pw = [0]*NORD; pw[0] = 1
    for c in coeffs:
        c = int(c) % P
        if c:
            for i in range(NORD):
                if pw[i]: out[i] = (out[i] + c*pw[i]) % P
        pw = smul(pw, tser)
    return out

Ns, Ds = polyser(N_QV), polyser(D_QV)
dNs, dDs = polyser(sp.diff(N_QV, t)), polyser(sp.diff(sp.expand(D_QV), t))
W1 = ddt(W); W2 = ddt(W1); W3 = ddt(W2)
lhs = smul(smul(Ds, Ds), W3)
tmp = smul(smul(Ns, Ds), W1)
lhs = [(lhs[i] + 4*tmp[i]) % P for i in range(NORD)]
tmp = [(smul(dNs, Ds)[i] - smul(Ns, dDs)[i]) % P for i in range(NORD)]
tmp = smul(tmp, W)
lhs = [(lhs[i] + 2*tmp[i]) % P for i in range(NORD)]
check("Ntilde(W) = 0 through q^%d (pole-free polynomial form)" % (NORD-7),
      all(x % P == 0 for x in lhs[:NORD-6]))

# ------------- CHECK 2: no P(t, t*W'/W) = 0 up to bidegree (28, 12)
Lam = divW(smul(tser, W1))                               # t * W'/W, unit series
check("Lambda = t*W'/W is a unit power series", Lam[0] % P == 1)

DT, DX, NEQ = 28, 12, 480
tp = [[0]*NORD for _ in range(DT+1)]; tp[0][0] = 1
for i in range(1, DT+1): tp[i] = smul(tp[i-1], tser)
Lp = [[0]*NORD for _ in range(DX+1)]; Lp[0][0] = 1
for j in range(1, DX+1): Lp[j] = smul(Lp[j-1], Lam)
cols = [smul(tp[i], Lp[j])[:NEQ] for i in range(DT+1) for j in range(DX+1)]

ncols = len(cols)
rows = [[cols[c][rr] for c in range(ncols)] for rr in range(NEQ)]
rank = 0
for col in range(ncols):
    piv = next((rr for rr in range(rank, NEQ) if rows[rr][col] % P), None)
    if piv is None:
        continue
    rows[rank], rows[piv] = rows[piv], rows[rank]
    inv = pow(rows[rank][col], P-2, P)
    rows[rank] = [(x*inv) % P for x in rows[rank]]
    for rr in range(NEQ):
        if rr != rank and rows[rr][col]:
            f = rows[rr][col]
            rows[rr] = [(a - f*b) % P for a, b in zip(rows[rr], rows[rank])]
    rank += 1
check("no relation P(t, tW'/W) = 0 with deg_t <= %d, deg_X <= %d" % (DT, DX),
      rank == ncols, "rank %d / %d unknowns, %d equations" % (rank, ncols, NEQ))

print()
if fails:
    print("RESULT: FAIL —", ", ".join(fails)); sys.exit(1)
print("RESULT: ALL CHECKS PASS")
sys.exit(0)
