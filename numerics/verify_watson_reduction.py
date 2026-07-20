"""
Fable's Watson-integral reduction of the hyperkagome/K4 LGF -- exact verification.

Chain:
 (1) D(t,kappa) := det(I4 - t * S_prim^2)  ==  (1 - 6t + (3-Lam)t^2)^2 - t^3 * Xi^2,
     with Lam = uv+vs+su+c.c., Xi = u+v+s+uvs+c.c.   [difference of squares, exact]
 (2) On the 4:1 covering torus  kappa = B theta,  B = [[1,-1,1],[1,1,-1],[-1,1,1]] (det 4):
     Lam = 2(cos2th1+cos2th2+cos2th3)   [by construction: k_i+k_j = 2 th_m]
     Xi  = 8 cos th1 cos th2 cos th3     [product identity -- verified symbolically]
 (3) Hence D = D+ * D-,  D+- = (1-6t+3t^2) - 2t^2*Sum cos2th_i -+ 8 t^{3/2} Prod cos th_i,
     and D-(theta) = D+(theta + (pi,pi,pi)):  CT log D- = CT log D+.
 (4) Phi(t) = (1/6) CT[4 - t d/dt log D] = 2/3 - (t/3) * d/dt CT_theta log D+.
     Series check against nu.json through m = 7.
==> the hyperkagome LGF is the LGF of the bcc lattice with 1st (weight t^{3/2}) and
    2nd (sc, weight t^2) neighbours at spectral parameter 1-6t+3t^2: a generalized
    Watson integral (Joyce class).
"""
import json, os, sys
from fractions import Fraction as F
import sympy as sp

u, v, s, t = sp.symbols('u v s t')

# ---------- (1) the difference-of-squares identity ----------
S = sp.Matrix([[0, 1, 1, 1],
               [1, 0, 1/v, u],
               [1, v, 0, 1/s],
               [1, 1/u, s, 0]])
D = sp.expand(sp.cancel(sp.together((sp.eye(4) - t * (S * S)).det())))
Lam = u*v + v*s + s*u + 1/(u*v) + 1/(v*s) + 1/(s*u)
Xi = u + v + s + u*v*s + 1/u + 1/v + 1/s + 1/(u*v*s)
Fform = (1 - 6*t + (3 - Lam)*t**2)**2 - t**3 * Xi**2
ok1 = sp.simplify(sp.together(D - sp.expand(Fform))) == 0
print("(1) det(I - t S^2) == (1-6t+(3-Lam)t^2)^2 - t^3 Xi^2 :", "PASS" if ok1 else "FAIL")

# Hermiticity sanity: S^dagger (u->1/u etc, transpose) equals S
Sd = S.T.subs({u: 1/u, v: 1/v, s: 1/s}, simultaneous=True)
print("    S Hermitian on the torus:", "PASS" if sp.simplify(Sd - S) == sp.zeros(4) else "FAIL")

# ---------- (2) the covering identities ----------
t1, t2, t3 = sp.symbols('t1 t2 t3', real=True)
k1, k2, k3 = t1 - t2 + t3, t1 + t2 - t3, -t1 + t2 + t3      # kappa = B theta ; k_i+k_j = 2 th
BB = sp.Matrix([[1, -1, 1], [1, 1, -1], [-1, 1, 1]])
print("(2) det B =", BB.det(), "(4:1 measure-preserving covering)")
subs_exp = {u: sp.exp(sp.I*k1), v: sp.exp(sp.I*k2), s: sp.exp(sp.I*k3)}
Lam_th = sp.simplify(sp.expand(sp.expand_trig(sp.re(sp.expand(Lam.subs(subs_exp).rewrite(sp.cos))))))
Xi_th = sp.simplify(sp.expand(sp.expand_trig(sp.re(sp.expand(Xi.subs(subs_exp).rewrite(sp.cos))))))
Lam_target = 2*(sp.cos(2*t1) + sp.cos(2*t2) + sp.cos(2*t3))
Xi_target = 8*sp.cos(t1)*sp.cos(t2)*sp.cos(t3)
ok2a = sp.simplify(sp.expand_trig(Lam_th - Lam_target)) == 0
ok2b = sp.simplify(sp.expand_trig(Xi_th - Xi_target)) == 0
print("    Lam(B theta) == 2 Sum cos 2th_i :", "PASS" if ok2a else "FAIL")
print("    Xi (B theta) == 8 Prod cos th_i :", "PASS" if ok2b else "FAIL")

# ---------- (3)+(4) series of the factored CT formula vs nu.json ----------
# Laurent-polynomial CT arithmetic in z_i = e^{i th_i}:  P1 = 2 Sum cos2th (sc),
# P2 = 8 Prod cos th (bcc).  D+ = g(t) - t^2 P1 - t^{3/2} P2,  g = 1-6t+3t^2.
# CT log D+ = log g + CT log(1 - Y/g),  Y = t^2 P1 + t^{3/2} P2  (h := sqrt(t)).
def lmul(a, b):
    r = {}
    for ka, ca in a.items():
        for kb, cb in b.items():
            k = (ka[0]+kb[0], ka[1]+kb[1], ka[2]+kb[2])
            r[k] = r.get(k, F(0)) + ca*cb
    return r
def lscale(a, c):
    return {k: cv*c for k, cv in a.items()}
def ladd(a, b):
    r = dict(a)
    for k, cv in b.items():
        r[k] = r.get(k, F(0)) + cv
        if r[k] == 0:
            del r[k]
    return r
P1 = {}
for i in range(3):
    for sgn in (2, -2):
        k = [0, 0, 0]; k[i] = sgn
        P1[tuple(k)] = F(1)                        # sum z_i^2 + z_i^-2  == 2 Sum cos2th
P2 = {}
for e1 in (1, -1):
    for e2 in (1, -1):
        for e3 in (1, -1):
            P2[(e1, e2, e3)] = F(1)                # Prod (z_i + 1/z_i) == 8 Prod cos th
MMAX = 8                                           # t-order to check (h-order 2*MMAX)
H = 2*MMAX + 1                                     # powers of h = sqrt(t): Y ~ h^3, h^4
# Y as dict: h-power -> Laurent poly
Ypow = {0: {(0, 0, 0): F(1)}}                      # Y^0 = 1
Y = {4: P1, 3: P2}
gser = [F(1)]                                      # 1/g series in t: g = 1-6t+3t^2
for n in range(1, MMAX+2):
    gser.append(6*gser[n-1] - 3*(gser[n-2] if n >= 2 else 0))
CTlog = {}                                         # h-power -> Fraction, CT log(1 - Y/g)
Yn = {0: {(0, 0, 0): F(1)}}
for n in range(1, H//3 + 2):
    newY = {}
    for hp, poly in Yn.items():
        for hy, py in Y.items():
            if hp + hy <= H:
                newY[hp+hy] = ladd(newY.get(hp+hy, {}), lmul(poly, py))
    Yn = newY                                      # Y^n
    # (1/g)^n as t-series: coefficients c with (1/g)^n = sum c_m t^m
    gn = [F(1)]
    for _ in range(n):
        gn = [sum(gn[j]*gser[m-j] for j in range(len(gn)) if 0 <= m-j < len(gser))
              for m in range(MMAX+1)]
    for hp, poly in Yn.items():
        ct = poly.get((0, 0, 0), F(0))
        if ct:
            for m, c in enumerate(gn):
                tot = hp + 2*m
                if tot <= H and c:
                    CTlog[tot] = CTlog.get(tot, F(0)) - F(1, n)*ct*c
half_ok = all(CTlog.get(p, F(0)) == 0 for p in range(1, H, 2))
print("(3) all half-integer t-powers of CT log D+ vanish:", "PASS" if half_ok else "FAIL")
# CT log D+ as t-series  = log g + collected even h-powers
loggser = [F(0)]                                   # log(1-6t+3t^2) = -sum (6t-3t^2)^n/n
for m in range(1, MMAX+1):
    loggser.append(F(0))
w = [F(0), F(6), F(-3)]
wn = [F(1)]
for n in range(1, MMAX+1):
    wn = [sum(wn[j]*w[m-j] for j in range(len(wn)) if 0 <= m-j < len(w))
          for m in range(MMAX+1)]
    for m, c in enumerate(wn):
        if m <= MMAX:
            loggser[m] -= F(1, n)*c if False else -F(1, n)*c*(-1)  # log(1-x): -x^n/n with x=6t-3t^2
# fix: log g = sum_n -(1/n) (6t-3t^2)^n  -> recompute cleanly
loggser = [F(0)]*(MMAX+1)
wn = [F(1)]
for n in range(1, MMAX+1):
    wn = [sum(wn[j]*w[m-j] for j in range(len(wn)) if 0 <= m-j < len(w))
          for m in range(MMAX+1)]
    for m, c in enumerate(wn):
        loggser[m] -= F(1, n)*c
J = [loggser[m] + CTlog.get(2*m, F(0)) for m in range(MMAX+1)]   # CT log D+ , t^m coeffs
# Phi = 2/3 - (t/3) d/dt J  ->  Phi_m = -(m/3) J_m  for m>=1 ; Phi_0 = 2/3
nu = [F(a, b) for a, b in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "nu.json")))]
pred = [F(2, 3)] + [-F(m, 3)*J[m] for m in range(1, MMAX+1)]
ok4 = all(pred[m] == nu[m] for m in range(MMAX+1))
print("(4) Phi(t) = 2/3 - (t/3) d/dt CT log D+  matches nu.json through m=%d :" % MMAX,
      "PASS" if ok4 else "FAIL")
print("    predicted:", [str(x) for x in pred])
print("    nu.json  :", [str(x) for x in nu[:MMAX+1]])

print()
if ok1 and ok2a and ok2b and half_ok and ok4:
    print("ALL PASS => the hyperkagome LGF is EXACTLY the generalized Watson integral of the")
    print("bcc lattice with 1st- (weight t^{3/2}) and 2nd- (sc, weight t^2) neighbour hopping")
    print("at spectral parameter 1-6t+3t^2:")
    print("   Phi(t) = 2/3 - (t/3) d/dt CT_theta log[ (1-6t+3t^2) - 2t^2 Sum_i cos 2th_i")
    print("                                            - 8 t^{3/2} cos th1 cos th2 cos th3 ]")
else:
    sys.exit(1)
