#!/usr/bin/env python3
"""
verify_moment_bridge.py  --  the missing link between the RAW moment series in
x = 1/E and the symmetry-reduced operator M in t = (E-1)^{-2}.

Context.  J.-M. Maillard (private communication, 2026-09-01) started from the raw
moment list in numerics/moments230.json,

    S(x) = sum_n m_n x^n = 1 + 4x^2 + 4x^3 + 28x^4 + 60x^5 + ... ,

guessed its minimal annihilator, and obtained an order-FIVE operator

    LCLM( N1 , NN3 . NN1 ),      N1  = Dx + (4x-1)/(2x^2-x-1),
                                 NN1 = Dx + 1/(x-1),
                                 NN3 = an order-three operator of degree 28,

and asked how NN3 relates to our L3 = M (order three, degree fifteen, in t).
Our README compresses the answer into eight words -- "after removing the
flat-band pole and using the exact reflection symmetry about E = 1" -- which is
exactly the clause he flagged as devilish.  This script makes the whole bridge
explicit and machine-checked, in exact rational arithmetic.

THE BRIDGE (everything below is verified, not asserted):

    z = E,   x = 1/z,   zeta = z - 1,   t = zeta^{-2} = x^2/(1-x)^2 .

    (A)  S(x) = (1/3)/(1+2x)  +  Phi( x^2/(1-x)^2 ) / (1-x)                [CHECK 2]

         first term  = the four flat bands at E = -2, spectral weight 1/3;
         second term = zeta*G_disp = Phi(t), regauged to the variable x.

    (B)  t = x^2/(1-x)^2 is the degree-two invariant of the Moebius INVOLUTION
         x -> x/(2x-1), i.e. of the reflection E -> 2-E.  So the passage
         S -> Phi is a quadratic pullback, not a substitution.               [CHECK 4]

    (C)  NN1 is NOT a delta function at E = 1.  It is the image of the trivial
         right factor d/dt of L4 = M.d/dt under the pullback + the gauge
         1/(1-x):  the constant solution of d/dt becomes 1/(1-x).            [CHECK 6]

    (D)  NN3 = (rational multiple of)  P(M) . [ (1-x)^4 / (2x) ],
         where P(M) is the pullback of M under t = x^2/(1-x)^2.  Verified
         coefficient by coefficient against Maillard's printed NN3.          [CHECK 7]

Run:  python verify_moment_bridge.py [path/to/numerics]

It needs only three files from the public repository
(github.com/bryannasr4-gif/hyperkagome-lgf): numerics/moments230.json,
numerics/nu.json and numerics/M_coeffs.json.  Give their directory as the first
argument, or set HYPERKAGOME_NUMERICS, or run the script from anywhere inside a
checkout.  SymPy is the only dependency.

The primitives are validated on operators of known structure before the real
computation, and the substantive checks each carry a NEGATIVE CONTROL: a
deliberately wrong variant that must FAIL.  Exit 0 means everything passed.
"""
import json, os, sys
from fractions import Fraction as Fr
from math import comb

import sympy as sp

DATA = ("moments230.json", "nu.json", "M_coeffs.json")


def find_numerics():
    """Locate the three repository data files, wherever this script was put."""
    cands = []
    if len(sys.argv) > 1:
        cands.append(sys.argv[1])
    if os.environ.get("HYPERKAGOME_NUMERICS"):
        cands.append(os.environ["HYPERKAGOME_NUMERICS"])
    for base in (os.path.dirname(os.path.abspath(__file__)), os.getcwd()):
        root = base
        for _ in range(5):
            cands += [root,
                      os.path.join(root, "numerics"),
                      os.path.join(root, "hyperkagome-lgf", "numerics")]
            root = os.path.dirname(root)
    for c in cands:
        if c and all(os.path.exists(os.path.join(c, f)) for f in DATA):
            return c
    raise SystemExit("could not find %s -- pass their directory as the first "
                     "argument, e.g.  python %s /path/to/numerics"
                     % (", ".join(DATA), os.path.basename(__file__)))


NUM = find_numerics()

x = sp.Symbol('x')
t = sp.Symbol('t')

FAIL = []
def check(name, cond, extra=""):
    print("  %-68s %s%s" % (name, "PASS" if cond else "*** FAIL ***",
                            ("   " + extra) if extra else ""))
    if not cond:
        FAIL.append(name)
    return cond

# ---------------------------------------------------------------------------
# Ore-algebra primitives over Q(x):  an operator is a list [a0, a1, ...] meaning
# sum_i a_i (d/dx)^i.
# ---------------------------------------------------------------------------

def op_compose(A, B):
    """(sum a_i D^i) . (sum b_j D^j), with D^i b = sum_k C(i,k) b^(k) D^(i-k)."""
    res = [sp.Integer(0)] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        if a == 0:
            continue
        for j, b in enumerate(B):
            if b == 0:
                continue
            for k in range(i + 1):
                res[i - k + j] += a * sp.binomial(i, k) * sp.diff(b, x, k)
    return [sp.cancel(c) for c in res]


def op_mul_right_function(A, f):
    """A . (multiplication by f)."""
    return op_compose(A, [f])


def op_apply_series(A, ser, N):
    """Apply A (coefficients rational in x) to a power series given as a list of
    Fractions, returning the first N coefficients of the image.  Used for the
    exact series checks; coefficients are cleared to polynomials first."""
    A = normalize_op(A)
    polys = [sp.Poly(a, x) for a in A]
    out = [Fr(0)] * N
    for i, P in enumerate(polys):
        cf = [Fr(int(c)) for c in reversed(P.all_coeffs())]   # cf[j] = coeff of x^j
        for j, c in enumerate(cf):
            if c == 0:
                continue
            # x^j * D^i ser  ->  coefficient of x^k is c * ff(k-j+i, i) * ser[k-j+i]
            for k in range(N):
                n = k - j + i
                if 0 <= n < len(ser):
                    ff = 1
                    for a in range(i):
                        ff *= (n - a)
                    out[k] += c * ff * ser[n]
    return out


def normalize_op(A):
    """Canonical form: clear denominators, remove the polynomial content and the
    integer content, fix the sign from the leading coefficient of the top-order
    term.  Two operators define the same solution space-with-multiplicity iff
    their canonical forms agree."""
    A = [sp.cancel(sp.together(a)) for a in A]
    while A and A[-1] == 0:
        A = A[:-1]
    den = sp.Integer(1)
    for a in A:
        den = sp.lcm(den, sp.denom(a))
    nums = [sp.expand(sp.cancel(a * den)) for a in A]
    g = sp.Integer(0)
    for n in nums:
        g = sp.gcd(g, n)
    nums = [sp.expand(sp.cancel(n / g)) for n in nums]
    polys = [sp.Poly(n, x) for n in nums]
    cont = sp.Integer(0)
    for P in polys:
        cont = sp.gcd(cont, P.content() if P.total_degree() >= 0 else 0)
    polys = [sp.Poly([c / cont for c in P.all_coeffs()], x) for P in polys]
    if polys[-1].LC() < 0:
        polys = [sp.Poly([-c for c in P.all_coeffs()], x) for P in polys]
    return [P.as_expr() for P in polys]


def pullback_gauge(L_t, phi, m, var=t):
    """Given L_t = [b_0..b_n] in the variable `var`, the map var = phi(x) and a
    gauge m(x), return the operator in x annihilating  y(x) = m(x) F(phi(x))
    for every solution F of L_t F = 0.

    Method: y^(j) = sum_i A[j][i] F^(i)(phi) with A lower triangular and
    A[j][j] = m phi'^j; invert A by forward substitution and substitute."""
    n = len(L_t) - 1
    dphi = sp.cancel(sp.diff(phi, x))
    rows = []
    a = [m] + [sp.Integer(0)] * n
    for j in range(n + 1):
        rows.append(list(a))
        na = [sp.Integer(0)] * (n + 1)
        for i in range(n + 1):
            if a[i] == 0:
                continue
            na[i] += sp.diff(a[i], x)
            if i + 1 <= n:
                na[i + 1] += a[i] * dphi
        a = [sp.cancel(v) for v in na]
    # B = A^{-1}, both lower triangular:  F^(i) = sum_j B[i][j] y^(j)
    B = [[sp.Integer(0)] * (n + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        B[i][i] = sp.cancel(1 / rows[i][i])
        for j in range(i - 1, -1, -1):
            s = sum(rows[i][k] * B[k][j] for k in range(j, i))
            B[i][j] = sp.cancel(-s / rows[i][i])
    Lphi = [sp.cancel(b.subs(var, phi)) for b in L_t]
    out = []
    for j in range(n + 1):
        out.append(sp.cancel(sum(Lphi[i] * B[i][j] for i in range(n + 1))))
    return out


def op_rightdivide(A, B):
    """Right Euclidean division A = Q.B + R in Q(x)<D>, B monic in D of order 1.
    Returns (Q, R)."""
    assert sp.cancel(B[-1] - 1) == 0, "divisor must be monic"
    A = [sp.cancel(a) for a in A]
    dB = len(B) - 1
    Q = [sp.Integer(0)] * (len(A) - dB)
    while len(A) - 1 >= dB and any(c != 0 for c in A):
        n = len(A) - 1
        while n >= 0 and A[n] == 0:
            n -= 1
        if n < dB:
            break
        lead = A[n]
        k = n - dB
        Q[k] = sp.cancel(Q[k] + lead)
        sub = op_compose([sp.Integer(0)] * k + [lead], B)
        A = [sp.cancel(A[i] - (sub[i] if i < len(sub) else 0)) for i in range(max(len(A), len(sub)))]
        while len(A) > 1 and A[-1] == 0:
            A = A[:-1]
    return Q, A

# ---------------------------------------------------------------------------
print(__doc__.split("Run:")[0].strip()[:0] or "", end="")
print("=" * 78)
print("verify_moment_bridge.py -- raw moment series in x  <->  operator M in t")
print("=" * 78)

# ---------------------------------------------------------------------------
print("\n[0] VALIDATION OF THE PRIMITIVES ON OPERATORS OF KNOWN STRUCTURE")
# 0a. composition against a hand computation:  D . x = x D + 1
c = op_compose([sp.Integer(0), sp.Integer(1)], [x])
check("(0a) D . x  ==  x D + 1", sp.cancel(c[0] - 1) == 0 and sp.cancel(c[1] - x) == 0)

# 0b. pullback of D_t^2 (solutions 1, t) under t = phi must annihilate 1 and phi
phi_test = x ** 3 - 2 * x
P = pullback_gauge([sp.Integer(0), sp.Integer(0), sp.Integer(1)], phi_test, sp.Integer(1))
def op_apply_expr(A, f):
    return sp.cancel(sum(A[i] * sp.diff(f, x, i) for i in range(len(A))))
ok = (op_apply_expr(P, sp.Integer(1)) == 0 and op_apply_expr(P, phi_test) == 0)
check("(0b) pullback of D_t^2 under t=x^3-2x kills 1 and phi(x)", ok)
check("(0b) NEGATIVE CONTROL: the same operator does NOT kill x^2",
      op_apply_expr(P, x ** 2) != 0)

# 0c. gauge: pullback with gauge m must annihilate m*1 and m*phi
mtest = 1 / (1 + x)
Pg = pullback_gauge([sp.Integer(0), sp.Integer(0), sp.Integer(1)], phi_test, mtest)
ok = (op_apply_expr(Pg, mtest) == 0 and op_apply_expr(Pg, sp.cancel(mtest * phi_test)) == 0)
check("(0c) gauged pullback kills m and m*phi", ok)
check("(0c) NEGATIVE CONTROL: it does NOT kill phi alone", op_apply_expr(Pg, phi_test) != 0)

# 0d. right division: (A.B) / B == A, and a non-multiple leaves a remainder
Atest = [x ** 2, sp.Integer(1) + x, sp.Integer(1)]
Btest = [1 / (x - 1), sp.Integer(1)]
Q, R = op_rightdivide(op_compose(Atest, Btest), Btest)
check("(0d) right division recovers the quotient exactly",
      all(sp.cancel(Q[i] - Atest[i]) == 0 for i in range(len(Atest))) and all(r == 0 for r in R))
Q2, R2 = op_rightdivide([sp.Integer(1), sp.Integer(0), sp.Integer(0), sp.Integer(1)], Btest)
check("(0d) NEGATIVE CONTROL: a non-multiple leaves a nonzero remainder",
      any(r != 0 for r in R2))

# ---------------------------------------------------------------------------
print("\n[1] THE DATA (read from the repository at run time)")
m = json.load(open(os.path.join(NUM, "moments230.json")))
nu = [Fr(a, b) for a, b in json.load(open(os.path.join(NUM, "nu.json")))]
Mc = json.load(open(os.path.join(NUM, "M_coeffs.json")))
print("    moments230.json : %d integers, m0..m10 = %s" % (len(m), m[:11]))
print("    nu.json         : %d rationals, nu0..nu4 = %s" % (len(nu), [str(v) for v in nu[:5]]))
print("    M_coeffs.json   : order %d, degrees %s" %
      (Mc["order"], [len(c) - 1 for c in Mc["coeffs"]]))

M_t = [sum(sp.Integer(int(cij)) * t ** j for j, cij in enumerate(col)) for col in Mc["coeffs"]]
check("(1a) M has order three and degree fifteen",
      len(M_t) == 4 and max(sp.Poly(c, t).degree() for c in M_t) == 15)

# ---------------------------------------------------------------------------
print("\n[2] THE MOMENT DECOMPOSITION   S(x) = (1/3)/(1+2x) + Phi(t(x))/(1-x)")
N = len(m)                      # 231 coefficients, x^0 .. x^230

def series_t_of_x(N):
    """t = x^2/(1-x)^2 = sum_{n>=2} (n-1) x^n."""
    return [Fr(0)] * 2 + [Fr(n - 1) for n in range(2, N)]

def compose_phi(nu, tser, N):
    """Phi(t(x)) truncated at x^N, by accumulating powers of t."""
    out = [Fr(0)] * N
    out[0] = nu[0]
    pw = [Fr(0)] * N
    pw[0] = Fr(1)
    for mm in range(1, N):
        if 2 * mm >= N:
            break
        if mm >= len(nu):
            raise RuntimeError("nu series too short for x^%d" % N)
        # pw <- pw * tser
        new = [Fr(0)] * N
        for i in range(N):
            if pw[i] == 0:
                continue
            for j in range(2, N - i):
                if tser[j] == 0:
                    continue
                new[i + j] += pw[i] * tser[j]
        pw = new
        for i in range(N):
            if pw[i]:
                out[i] += nu[mm] * pw[i]
    return out

def mul_geom(ser, N):
    """multiply by 1/(1-x): partial sums."""
    out = [Fr(0)] * N
    acc = Fr(0)
    for i in range(N):
        acc += ser[i]
        out[i] = acc
    return out

tser = series_t_of_x(N)
PhiC = compose_phi(nu, tser, N)
disp = mul_geom(PhiC, N)
flat = [Fr(1, 3) * Fr((-2) ** n) for n in range(N)]
rebuilt = [flat[i] + disp[i] for i in range(N)]
agree = sum(1 for i in range(N) if rebuilt[i] == Fr(m[i]))
check("(2a) S(x) reproduced EXACTLY from Phi for x^0..x^%d (%d/%d)" % (N - 1, agree, N),
      agree == N)

# negative controls
bad1 = [Fr(1, 2) * Fr((-2) ** n) + disp[i] for i, n in enumerate(range(N))]
check("(2b) NEGATIVE CONTROL: flat-band weight 1/2 instead of 1/3 fails",
      sum(1 for i in range(N) if bad1[i] == Fr(m[i])) < N)
tser_bad = [Fr(0)] * 2 + [Fr((-1) ** n * (n - 1)) for n in range(2, N)]   # x^2/(1+x)^2
PhiCb = compose_phi(nu, tser_bad, N)
bad2 = [flat[i] + v for i, v in enumerate(mul_geom(PhiCb, N))]
check("(2c) NEGATIVE CONTROL: t = x^2/(1+x)^2 fails",
      sum(1 for i in range(N) if bad2[i] == Fr(m[i])) < N)
bad3 = [flat[i] + PhiC[i] for i in range(N)]
check("(2d) NEGATIVE CONTROL: omitting the gauge 1/(1-x) fails",
      sum(1 for i in range(N) if bad3[i] == Fr(m[i])) < N)

# ---------------------------------------------------------------------------
print("\n[3] THE SYMMETRY-REDUCED SERIES nu REBUILT FROM THE RAW MOMENTS")
d = [Fr(m[n]) - Fr(1, 3) * Fr((-2) ** n) for n in range(N)]
def central(p):
    return sum(Fr(comb(p, k)) * Fr((-1) ** (p - k)) * d[k] for k in range(p + 1))
even_ok = all(central(2 * mm) == nu[mm] for mm in range(len(nu)))
check("(3a) nu_m = int rho_disp (E-1)^{2m} matches nu.json for m=0..%d" % (len(nu) - 1), even_ok)
odd_zero = [p for p in range(1, N, 2) if central(p) != 0]
check("(3b) ALL odd central moments vanish, p = 1,3,...,%d (the E -> 2-E symmetry)" % (N - 2),
      odd_zero == [], "first nonzero: %s" % (odd_zero[:1] or "none"))
# negative control: the same test on the RAW moments (flat band not removed) must fail
draw = [Fr(m[n]) for n in range(N)]
def central_raw(p):
    return sum(Fr(comb(p, k)) * Fr((-1) ** (p - k)) * draw[k] for k in range(p + 1))
check("(3c) NEGATIVE CONTROL: without removing the flat band, odd moments do NOT vanish",
      any(central_raw(p) != 0 for p in range(1, 12, 2)))

# ---------------------------------------------------------------------------
print("\n[4] t = x^2/(1-x)^2 IS THE INVARIANT OF THE INVOLUTION E -> 2-E")
sigma = x / (2 * x - 1)                       # E -> 2-E  read in x = 1/E
check("(4a) x -> x/(2x-1) is an involution", sp.cancel(sigma.subs(x, sigma) - x) == 0)
check("(4b) it is the reflection E -> 2-E", sp.cancel(1 / sigma.subs(x, 1 / t) - (2 - t)) == 0)
phi = x ** 2 / (1 - x) ** 2
check("(4c) t = x^2/(1-x)^2 is sigma-invariant", sp.cancel(phi.subs(x, sigma) - phi) == 0)
check("(4d) the map has degree two", sp.degree(sp.numer(sp.cancel(phi)), x) == 2
      and sp.degree(sp.denom(sp.cancel(phi)), x) == 2)
# branch dictionary
print("    branch dictionary (t <- the two E-preimages):")
for tv, lab in [(sp.Rational(1, 9), "E = 4 (band edge) and E = -2 (flat band)"),
                (sp.Rational(1, 5), "E = 1 +- sqrt5"),
                (sp.Rational(1, 4), "E = 3 and E = -1"),
                (sp.Integer(1), "E = 2 and E = 0 (principal van Hove)"),
                (sp.Integer(0), "E = infinity (double)")]:
    sols = sp.solve(sp.Eq(phi, tv), x)
    print("      t = %-4s <-  x in %-28s   %s" % (tv, str(sols), lab))

# ---------------------------------------------------------------------------
print("\n[5] PULLING M BACK:  NN3 = P(M) . [(1-x)^4/(2x)]")
gauge = 1 / (1 - x)                      # S_disp(x) = Phi(t(x)) * gauge
P_M = pullback_gauge(M_t, phi, sp.Integer(1))          # annihilates Phi(t(x))
f = (1 - x) ** 4 / (2 * x)
NN3_ours = normalize_op(op_mul_right_function(P_M, f))
degs = [sp.Poly(c, x).degree() for c in NN3_ours]
print("    our NN3 : order %d, coefficient degrees %s" % (len(NN3_ours) - 1, degs))

L4_t = [sp.Integer(0)] + M_t                            # L4 = M . d/dt
P_L4 = pullback_gauge(L4_t, phi, gauge)                 # annihilates Phi(t(x))/(1-x)
N4 = normalize_op(P_L4)
print("    our N4  : order %d, coefficient degrees %s" %
      (len(N4) - 1, [sp.Poly(c, x).degree() for c in N4]))

NN1 = [1 / (x - 1), sp.Integer(1)]
Q, R = op_rightdivide(N4, NN1)
check("(5a) N4 is EXACTLY right-divisible by NN1 = Dx + 1/(x-1)", all(r == 0 for r in R))
check("(5b) N4 / NN1 == P(M).[(1-x)^4/(2x)]  (canonical forms agree)",
      normalize_op(Q) == NN3_ours)

# ---------------------------------------------------------------------------
print("\n[6] WHY NN1 APPEARS: it is the gauged trivial right factor d/dt of L4")
Dt = [sp.Integer(0), sp.Integer(1)]
P_Dt = pullback_gauge(Dt, phi, gauge)     # annihilates (const)*gauge
check("(6a) pullback+gauge of d/dt is a rational multiple of NN1",
      normalize_op(P_Dt) == normalize_op(NN1))
check("(6b) so the solution 1/(1-x) of NN1 is the gauged CONSTANT, not a delta at E=1",
      sp.cancel(op_apply_expr(NN1, gauge)) == 0)
# negative control: a delta at E = 1 would be an extra 1/(1-x) in S; show S has none
check("(6c) NEGATIVE CONTROL: S(x) needs no independent 1/(1-x) term "
      "(check 2a already closed the books)", agree == N)

# ---------------------------------------------------------------------------
print("\n[7] COMPARISON WITH MAILLARD'S PRINTED NN3 (order three, degree 28)")
P14 = (25920*x**14 - 67584*x**13 + 54880*x**12 - 488032*x**11 + 217416*x**10
       + 366368*x**9 - 110648*x**8 - 154392*x**7 - 1365*x**6 + 92072*x**5
       - 37156*x**4 - 1696*x**3 + 4024*x**2 - 896*x + 64)
P22 = (54743040*x**22 - 157994496*x**21 + 92272896*x**20 - 1159898880*x**19
       + 542944384*x**18 + 2539682944*x**17 - 1230378464*x**16 - 2079840288*x**15
       + 766774496*x**14 + 1278571312*x**13 - 375569060*x**12 - 639112648*x**11
       + 295033624*x**10 + 132779069*x**9 - 115662716*x**8 + 10603561*x**7
       + 12818606*x**6 - 4351956*x**5 - 1280*x**4 + 276040*x**3 - 65520*x**2
       + 6592*x - 256)
P24 = (296110080*x**24 - 1361055744*x**23 + 2386980864*x**22 - 10557751296*x**21
       + 18033853952*x**20 + 2804154368*x**19 - 22727726080*x**18 + 2408041728*x**17
       + 14110613504*x**16 + 3494831456*x**15 - 14226801152*x**14 + 1544149824*x**13
       + 7786156932*x**12 - 4284833964*x**11 - 489878130*x**10 + 1159823544*x**9
       - 337924079*x**8 - 36869110*x**7 + 43130027*x**6 - 8786544*x**5
       - 135380*x**4 + 336304*x**3 - 55624*x**2 + 3584*x - 64)
P25 = (348364800*x**25 - 1860433920*x**24 + 4225941504*x**23 - 17915424768*x**22
       + 35409525248*x**21 - 16308632576*x**20 - 8359990784*x**19 - 15977177344*x**18
       + 28290102336*x**17 + 17457912000*x**16 - 47218896448*x**15 + 19688934784*x**14
       + 13108308856*x**13 - 15560604020*x**12 + 4599412624*x**11 + 967410190*x**10
       - 985790086*x**9 + 162055029*x**8 + 73091664*x**7 - 44742005*x**6
       + 11267048*x**5 - 1689844*x**4 + 178080*x**3 - 16712*x**2 + 1408*x - 64)
NN3_jmm = [
    P25,
    P24 * x * (2*x - 1),
    x**2 * (x - 1) * (2*x - 1)**2 * P22,
    x**3 * (3*x - 1) * (2*x + 1) * (4*x - 1) * (x + 1) * (4*x**2 + 2*x - 1)
        * P14 * (x - 1)**2 * (2*x - 1)**3,
]
NN3_jmm_n = normalize_op(NN3_jmm)
check("(7a) OUR pulled-back NN3 == MAILLARD'S NN3, coefficient by coefficient",
      NN3_jmm_n == NN3_ours)
if NN3_jmm_n != NN3_ours:
    for i in range(4):
        print("       D^%d : %s" % (i, "equal" if sp.expand(NN3_jmm_n[i] - NN3_ours[i]) == 0
                                    else "DIFFER"))
# negative control: perturb one coefficient
bad = list(NN3_jmm)
bad[0] = bad[0] + 1
check("(7b) NEGATIVE CONTROL: perturbing one coefficient breaks the match",
      normalize_op(bad) != NN3_ours)

# the degree-14 factor is exactly the pullback of p_7
p7 = (101025*t**7 - 369600*t**6 + 455798*t**5 - 290956*t**4 + 93657*t**3
      - 17580*t**2 + 1800*t - 64)
p7_pull = sp.expand(sp.cancel(p7.subs(t, phi) * (1 - x) ** 14))
check("(7c) Maillard's degree-14 factor == -(1-x)^14 p_7(x^2/(1-x)^2)",
      sp.expand(p7_pull + P14) == 0)

# the head coefficient dictionary, factor by factor
head_pairs = [("(9t-1)", 9*t - 1, (4*x - 1) * (2*x + 1)),
              ("(5t-1)", 5*t - 1, 4*x**2 + 2*x - 1),
              ("(4t-1)", 4*t - 1, (3*x - 1) * (x + 1)),
              ("(t-1)", t - 1, 2*x - 1)]
allok = True
for lab, ex, expect in head_pairs:
    got = sp.expand(sp.cancel(ex.subs(t, phi) * (1 - x) ** 2))
    allok = allok and sp.expand(got - expect) == 0
check("(7d) each singular factor of M pulls back to the printed factor of NN3", allok)

# ---------------------------------------------------------------------------
print("\n[8] SERIES CHECK OF THE OPERATORS AGAINST THE RAW MOMENTS")
NSER = 200
Sd = disp[:NSER + 40]                             # the dispersive part of S(x)
img = op_apply_series(op_compose(NN3_jmm, NN1), Sd, NSER)
check("(8a) Maillard's NN3.NN1 annihilates Phi(t(x))/(1-x) through x^%d" % (NSER - 1),
      all(v == 0 for v in img))
img2 = op_apply_series(NN3_jmm, Sd, NSER)
check("(8b) NEGATIVE CONTROL: NN3 alone does NOT annihilate it",
      any(v != 0 for v in img2))
# N1 and the flat band
N1op = [(4*x - 1) / (2*x**2 - x - 1), sp.Integer(1)]
sol_N1 = 1 / ((2*x + 1) * (x - 1))
check("(8c) N1 annihilates 1/((2x+1)(x-1))", sp.cancel(op_apply_expr(N1op, sol_N1)) == 0)
check("(8d) 1/((2x+1)(x-1)) = -(2/3)/(2x+1) - (1/3)/(1-x): so modulo the NN1 solution "
      "1/(1-x)\n       it is the flat-band pole -- the order-one factor carries no extra "
      "physics",
      sp.cancel(sol_N1 - (sp.Rational(-2, 3) / (2*x + 1) + sp.Rational(-1, 3) / (1 - x))) == 0)


# ---------------------------------------------------------------------------
print("\n[9] CLOSING MAILLARD'S LOOP: the order-five LCLM annihilates the RAW series")
# LCLM(N1, A) = (D - v'/v) . A  with v = A[u], u the solution of N1.
A5 = op_compose(NN3_jmm, NN1)
u = sol_N1
v = sp.cancel(op_apply_expr(A5, u))
check("(9a) v = (NN3.NN1)[1/((2x+1)(x-1))] is a nonzero RATIONAL function",
      v != 0 and v.is_rational_function(x))
L5 = op_compose([sp.cancel(-sp.diff(v, x) / v), sp.Integer(1)], A5)
NS = 160
Sraw = [Fr(mm) for mm in m]
img5 = op_apply_series(L5, Sraw, NS)
check("(9b) LCLM(N1, NN3.NN1) annihilates S(x) = sum m_n x^n through x^%d" % (NS - 1),
      all(w == 0 for w in img5))
img4 = op_apply_series(A5, Sraw, NS)
check("(9c) NEGATIVE CONTROL: NN3.NN1 alone does NOT annihilate S(x) "
      "(his own observation, reproduced)", any(w != 0 for w in img4))
check("(9d) 1/(1+2x) = -(3/2)u - (1/2)/(1-x): Maple's N1 is the flat-band pole modulo\n"
      "       the gauged constant, i.e. a basis choice, not extra physics",
      sp.cancel(1 / (2*x + 1) - (sp.Rational(-3, 2) * u + sp.Rational(-1, 2) / (1 - x))) == 0)

# ---------------------------------------------------------------------------
print("\n[10] WHY NN3 IS HOMOMORPHIC TO ITS ADJOINT (his observation, explained)")
# Sym^2(M) has the rational solution R (numerics/CERTIFICATE_orthogonal.txt).  Under
# Y -> Y(phi(x))/f(x) every quadratic invariant transforms as R -> R(phi)/f^2, so:
R_t = (-sp.Rational(1, 272) * (15*t**2 + 17*t - 8)**2
       / (t**2 * (t - 1)**2 * (4*t - 1) * (5*t - 1) * (9*t - 1)))
R_x = sp.cancel(R_t.subs(t, phi) / f**2)
print("    predicted invariant of NN3:  R_NN3(x) = R(x^2/(1-x)^2) * 4x^2/(1-x)^8")


def series_of_ratfun(expr, x0, N):
    X = sp.Symbol('X')
    e = sp.cancel(sp.together(expr.subs(x, x0 + X)))
    nu_, de_ = sp.fraction(e)
    num = [Fr(int(sp.Rational(c).p), int(sp.Rational(c).q))
           for c in sp.Poly(sp.expand(nu_), X).all_coeffs()[::-1]] + [Fr(0)] * N
    den = [Fr(int(sp.Rational(c).p), int(sp.Rational(c).q))
           for c in sp.Poly(sp.expand(de_), X).all_coeffs()[::-1]] + [Fr(0)] * N
    assert den[0] != 0
    out = [Fr(0)] * N
    for n in range(N):
        out[n] = (num[n] - sum(den[k] * out[n - k] for k in range(1, n + 1))) / den[0]
    return out


def local_solutions(A, x0, N):
    """The three Taylor solutions of an order-3 operator at an ordinary point x0."""
    X = sp.Symbol('X')
    a = []
    for ai in A:
        P = sp.Poly(sp.expand(ai.subs(x, x0 + X)), X)
        a.append([Fr(int(sp.Rational(c).p), int(sp.Rational(c).q))
                  for c in P.all_coeffs()[::-1]])
    lead = a[3][0]
    assert lead != 0, "x0 is a singular point"
    sols = []
    for ic in range(3):
        c = [Fr(0)] * N
        c[ic] = Fr(1)
        for k in range(0, N - 3):
            tot = Fr(0)
            for i in range(4):
                for j, aij in enumerate(a[i]):
                    if aij == 0 or (i == 3 and j == 0):
                        continue
                    n = k - j + i
                    if 0 <= n < N:
                        ff = 1
                        for q in range(i):
                            ff *= (n - q)
                        tot += aij * ff * c[n]
            c[k + 3] = -tot / (lead * ((k + 3) * (k + 2) * (k + 1)))
        sols.append(c)
    return sols


x0 = sp.Rational(1, 10)
check("(10a) x0 = 1/10 is an ordinary point of NN3",
      sp.cancel(NN3_jmm_n[3].subs(x, x0)) != 0)
NLOC = 40
ys = local_solutions(NN3_jmm_n, x0, NLOC)
NN3_shift = [sp.expand(ai.subs(x, x0 + x)) for ai in NN3_jmm_n]
check("(10a') the local solver really solves NN3 (residual of y_1 vanishes)",
      all(w == 0 for w in op_apply_series(NN3_shift, ys[0], NLOC - 6)))
check("(10a'') NEGATIVE CONTROL: a perturbed series is not a solution",
      any(w != 0 for w in op_apply_series(NN3_shift, [c + Fr(1) for c in ys[0]], NLOC - 6)))
Rser = series_of_ratfun(R_x, x0, NLOC)
pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
rows = []
for k in range(NLOC):
    row = []
    for (i, j) in pairs:
        val = sum(ys[i][p] * ys[j][k - p] for p in range(k + 1))
        row.append(val * (1 if i == j else 2))
    rows.append(row)


def fit_and_count(target):
    Amat = sp.Matrix([[sp.Rational(vv.numerator, vv.denominator) for vv in r] for r in rows[:6]])
    bvec = sp.Matrix([sp.Rational(vv.numerator, vv.denominator) for vv in target[:6]])
    cs = Amat.solve(bvec)
    bad = 0
    for k in range(NLOC):
        val = sum(sp.Rational(rows[k][q].numerator, rows[k][q].denominator) * cs[q]
                  for q in range(6))
        if sp.nsimplify(val - sp.Rational(target[k].numerator, target[k].denominator)) != 0:
            bad += 1
    return bad, cs


nbad, cs = fit_and_count(Rser)
check("(10b) R_NN3 = sum c_ij y_i y_j with CONSTANT c_ij: %d relations vs 6 unknowns, "
      "margin %d" % (NLOC, NLOC - 6), nbad == 0)
nbad_b, _ = fit_and_count(series_of_ratfun(sp.cancel(R_x * (1 + x)), x0, NLOC))
check("(10c) NEGATIVE CONTROL: R_NN3*(1+x) is NOT such a quadratic invariant", nbad_b != 0)
print("    => Sym^2(NN3) has a rational solution, so NN3 preserves a nondegenerate symmetric")
print("       form and is homomorphic to its adjoint -- exactly what he observed.  It is")
print("       inherited from M: pullback and gauge send R(t) to R(phi(x))/f(x)^2.")

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
if FAIL:
    print("*** %d CHECK(S) FAILED: %s" % (len(FAIL), FAIL))
    sys.exit(1)
print("ALL CHECKS PASS.")
print("""
SUMMARY OF THE BRIDGE

  x = 1/E,  zeta = E-1,  t = zeta^-2 = x^2/(1-x)^2 .

  S(x) = sum m_n x^n = (1/3)/(1+2x)  +  Phi(t(x))/(1-x)

  * (1/3)/(1+2x) is the flat band at E = -2, weight (12-8)/12 = 1/3.  It is the
    piece removed by d_n = m_n - (1/3)(-2)^n.
  * t = x^2/(1-x)^2 is the degree-two invariant of the involution x -> x/(2x-1),
    which is the reflection E -> 2-E.  So the reduction is a QUADRATIC PULLBACK,
    which is why one order-three operator in t produces an order-three operator
    of twice the degree in x, with each singular point of M splitting into its
    two E-preimages.
  * the gauge 1/(1-x) is zeta/x = the Jacobian of writing zeta*G_disp instead
    of G_disp.  Its constant carries the solution 1/(1-x): that is NN1, the
    pullback of the trivial right factor d/dt of L4 = M.d/dt.  There is NO
    Dirac mass at E = 1.

  NN3 = (rational multiple of)  P(M) . [(1-x)^4/(2x)],  P = pullback along t(x).
""")
