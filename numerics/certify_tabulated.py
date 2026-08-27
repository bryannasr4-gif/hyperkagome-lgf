"""
certify_tabulated.py  --  exact certificate that the uniformizing equation of V2 is the entry
                          already printed in two published tables of genus-zero groups.

The projective normal form of V2 is  w'' + Q_V(t) w = 0  with

  Q_V = N(t) / [4 t^2 (t-1)^2 (4t-1)^2 (5t-1)^2 (9t-1)^2],
  N   = 24300 t^8 - 58860 t^7 + 73437 t^6 - 44294 t^5 + 15111 t^4 - 3160 t^3 + 407 t^2 - 30 t + 1,

read here from numerics/V2_data.json (the certified V2 of the repository), not retyped.

Two published tables carry this object in other generators.  A Schwarzian potential is a
QUADRATIC DIFFERENTIAL, not a function: under a change of generator it transports with the
square of the Jacobian.  Comparing the tables by naive substitution therefore produces a
spurious mismatch, and both directions are certified below.

  PART 1 -- Lian and Yau, "Mirror maps, modular relations and hypergeometric series I",
            arXiv:hep-th/9507151, table of genus-zero groups, row Gamma_0(30)+.  The row
            prints an order-three operator L (in the Euler operator Theta = x d/dx) and a
            potential Q(x); both are transcribed here from the source.

            [1a] the printed Theta-form operator, converted to d/dx form, is the order-three
                 operator L3 used below -- coefficient by coefficient, no scalar factor;
            [1b] L3 = Sym^2(L2) for an explicit order-two L2;
            [1d] the printed Q(x) is the projective normal form of L2;
            [1e] pulled back by the Moebius map X = x/(1-x) as a quadratic differential,
                 the printed Q(x) is EXACTLY Q_V;
            [1h] the singular sets correspond, {0,1/8,1/4,1/3,-1,oo} -> {0,1/9,1/5,1/4,oo,1}.

  PART 2 -- Lian and Wiczer, "Genus zero modular functions", arXiv:math/0611291: the
            Schwarzian Q-value tabulated for Conway-Norton class 30B, in the generator
            z = 1/T_30B related to ours by 1/t = T_30B + 4.

            [2b] Q_V(t) dt^2 = Q_30B(z) dz^2 identically, with t = z/(1+4z);
            [2d] the naive substitution Q_V(t(z)) differs from Q_30B(z) by exactly (1+4z)^4;
            [2i] pole dictionary t = 0,1/9,1/5,1/4,1,oo <-> z = 0,1/5,1,oo,-1/3,-1/4;
            [2j] identical Laurent heads: 1/4 at the cusp, 3/16 at each order-two elliptic image.

  PART 3 -- consistency: the two tables are related to each other by the composed Moebius map
            X = z/(1+3z), computed here from the two individual maps.

The Sym^2 and Theta-to-D primitives are validated on operators of known structure BEFORE they
are used on the tabulated row, and every structural check carries a negative control that is
shown to fail.  What the tables do not carry is a proof: they are obtained by fitting an ansatz
to a bounded q-expansion.  The a priori pole-degree bound in the proof of the modularity theorem
is what converts such a match into an identity, and hence proves the tabulated entries.

Run:  python certify_tabulated.py > CERTIFICATE_tabulated.txt      -> exit 0 on success.
"""
import json
import os
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))

x, z, s = sp.symbols('x z s')

FAILURES = []


def check(name, ok):
    print("  [%s] %s" % ("PASS" if ok else "FAIL", name))
    if not ok:
        FAILURES.append(name)


def sym2(p, q, var):
    """Sym^2(D^2 + p D + q) = D^3 + 3p D^2 + (2p^2 + p' + 4q) D + (4pq + 2q')."""
    return (3*p,
            2*p**2 + sp.diff(p, var) + 4*q,
            4*p*q + 2*sp.diff(q, var))


def theta_to_D(theta_coeffs_by_power, var):
    """Convert  sum_i var^i * P_i(Theta),  Theta = var d/dvar,  to d/dvar form.

    theta_coeffs_by_power maps i -> list of coefficients of P_i in ascending powers of Theta.
    Returns the four coefficients of D^3, D^2, D^1, D^0 as a dict.
    """
    f = sp.Function('f')(var)
    expr = sp.Integer(0)
    for i, poly in theta_coeffs_by_power.items():
        term, cur = sp.Integer(0), f
        for k, c in enumerate(poly):
            if k:
                cur = sp.expand(var*sp.diff(cur, var))
            term += c*cur
        expr += var**i*term
    expr = sp.expand(expr)
    out = {k: sp.factor(sp.expand(expr.coeff(sp.Derivative(f, (var, k))))) for k in (3, 2, 1)}
    out[0] = sp.factor(sp.expand(expr.coeff(f)))
    return out


def tpoly(*factors):
    """Product of polynomials in Theta, each given in ascending-coefficient form."""
    r = [sp.Integer(1)]
    for g in factors:
        out = [sp.Integer(0)]*(len(r) + len(g) - 1)
        for i, ri in enumerate(r):
            for j, gj in enumerate(g):
                out[i + j] += ri*gj
        r = out
    return r


def head(expr, var, p, order=2):
    """Laurent head at a finite point: the coefficient of (var-p)^-order."""
    e = sp.Symbol('e')
    return sp.simplify(sp.limit(sp.together(expr.subs(var, p + e)*e**order), e, 0))


print("=" * 78)
print("certify_tabulated.py  --  V2 is the tabulated genus-zero entry Gamma_0(30)+")
print("=" * 78)
print()

# ---------------------------------------------------------------- Q_V, read from V2_data.json
with open(os.path.join(HERE, "V2_data.json")) as fh:
    V2 = json.load(fh)
QV = sp.simplify(sp.sympify(V2["Q"]["num"])/sp.sympify(V2["Q"]["den"])).subs(sp.Symbol('t'), x)

N_paper = (24300*x**8 - 58860*x**7 + 73437*x**6 - 44294*x**5
           + 15111*x**4 - 3160*x**3 + 407*x**2 - 30*x + 1)
QV_paper = N_paper/(4*x**2*(x - 1)**2*(4*x - 1)**2*(5*x - 1)**2*(9*x - 1)**2)

print("Check 0: primitives validated on objects of known structure")
check("Q_V read from V2_data.json equals equation (QV) as printed in the paper",
      sp.simplify(sp.together(QV - QV_paper)) == 0)

# Sym^2 on a control: p = 0 must give the form D^3 + 4q D + 2q' printed in the paper.
qc = sp.Function('qc')(x)
c2, c1, c0 = sym2(sp.Integer(0), qc, x)
check("CONTROL Sym^2(D^2 + q) == D^3 + 4q D + 2q'",
      sp.simplify(c2) == 0 and sp.simplify(c1 - 4*qc) == 0
      and sp.simplify(c0 - 2*sp.diff(qc, x)) == 0)

# Sym^2 of the Euler operator D^2 - (1/x)D (solutions 1, x^2) must kill the products {1,x^2,x^4}.
e2, e1, e0 = sym2(-1/x, sp.Integer(0), x)
ok = all(sp.simplify(sp.diff(m, x, 3) + e2*sp.diff(m, x, 2) + e1*sp.diff(m, x) + e0*m) == 0
         for m in (sp.Integer(1), x**2, x**4))
check("CONTROL Sym^2(D^2 - (1/x)D) annihilates {1, x^2, x^4}", ok)
check("NEGATIVE CONTROL that same operator does NOT annihilate x^3",
      sp.simplify(sp.diff(x**3, x, 3) + e2*sp.diff(x**3, x, 2)
                  + e1*sp.diff(x**3, x) + e0*x**3) != 0)

# Theta-to-D converter on a control: Theta^3 = x^3 D^3 + 3x^2 D^2 + x D.
T3 = theta_to_D({0: [0, 0, 0, 1]}, x)
check("CONTROL Theta^3 == x^3 D^3 + 3x^2 D^2 + x D",
      sp.expand(T3[3] - x**3) == 0 and sp.expand(T3[2] - 3*x**2) == 0
      and sp.expand(T3[1] - x) == 0 and sp.expand(T3[0]) == 0)
check("NEGATIVE CONTROL Theta^3 is NOT x^3 D^3 (the lower-order terms are real)",
      sp.expand(T3[2]) != 0)
print()

# ================================================================ PART 1: the Lian-Yau row
print("Check 1: Lian-Yau (arXiv:hep-th/9507151), table of genus-zero groups, row Gamma_0(30)+")

# The printed operator, transcribed from the source in the Euler operator Theta = x d/dx:
#   Theta^3 - 14 x^3 (1+T)(2+T)(3+2T) - 24 x^4 (2+T)(3+2T)(5+2T)
#           -      x (1+2T)(4+7T+7T^2) +    x^2 (1+T)(72+106T+53T^2)
one, two, three, five = (sp.Integer(k) for k in (1, 2, 3, 5))
printed = theta_to_D({
    0: [0, 0, 0, one],
    3: [-14*c for c in tpoly([one, one], [two, one], [three, two])],
    4: [-24*c for c in tpoly([two, one], [three, two], [five, two])],
    1: [-c for c in tpoly([one, two], [sp.Integer(4), sp.Integer(7), sp.Integer(7)])],
    2: [c for c in tpoly([one, one], [sp.Integer(72), sp.Integer(106), sp.Integer(53)])],
}, x)

# The same operator written directly in d/dx form.
L3 = {3: -x**3*(4*x - 1)*(3*x - 1)*(8*x - 1)*(x + 1),
      2: -3*x**2*(288*x**4 + 70*x**3 - 106*x**2 + 21*x - 1),
      1: -x*(1800*x**4 + 336*x**3 - 390*x**2 + 50*x - 1),
      0: -4*x*(180*x**3 + 21*x**2 - 18*x + 1)}
for k in (3, 2, 1, 0):
    check("[1a] printed Theta-form operator equals L3 in the D^%d coefficient" % k,
          sp.expand(printed[k] - L3[k]) == 0)
print("       leading coefficient: %s" % sp.factor(L3[3]))
check("NEGATIVE CONTROL a 1-unit change to the printed x^2 block breaks [1a]",
      sp.expand(theta_to_D({0: [0, 0, 0, one],
                            3: [-14*c for c in tpoly([one, one], [two, one], [three, two])],
                            4: [-24*c for c in tpoly([two, one], [three, two], [five, two])],
                            1: [-c for c in tpoly([one, two],
                                                  [sp.Integer(4), sp.Integer(7), sp.Integer(7)])],
                            2: [c for c in tpoly([one, one],
                                                 [sp.Integer(73), sp.Integer(106), sp.Integer(53)])],
                            }, x)[0] - L3[0]) != 0)

# The order-two operator whose symmetric square that row is.
den_LY = x*(x + 1)*(3*x - 1)*(4*x - 1)*(8*x - 1)
p_LY = (288*x**4 + 70*x**3 - 106*x**2 + 21*x - 1)/den_LY
q_LY = 2*(45*x**3 + 7*x**2 - 9*x + 1)/den_LY

s2, s1, s0 = sym2(p_LY, q_LY, x)
for k, (got, want) in enumerate([(s0, L3[0]/L3[3]), (s1, L3[1]/L3[3]), (s2, L3[2]/L3[3])]):
    check("[1b] Sym^2(L2) matches the tabulated operator in the D^%d coefficient" % k,
          sp.simplify(sp.together(got - want)) == 0)
check("NEGATIVE CONTROL a 1-unit perturbation of q_LY breaks the Sym^2 match",
      sp.simplify(sp.together(sym2(p_LY, q_LY + 1, x)[1] - L3[1]/L3[3])) != 0)

# The printed potential.
Q_table = ((1 - 22*x + 225*x**2 - 1292*x**3 + 4436*x**4 - 8304*x**5
            + 5124*x**6 + 2016*x**7 + 6912*x**8)
           / (4*(1 - 8*x)**2*(1 - 4*x)**2*(-1 - x)**2*x**2*(-1 + 3*x)**2))
Q_LY = sp.simplify(q_LY - p_LY**2/4 - sp.diff(p_LY, x)/2)
check("[1d] the printed potential Q(x) is the projective normal form of L2",
      sp.simplify(sp.together(Q_table - Q_LY)) == 0)

g = x/(1 - x)                                    # the Moebius change of generator X = x/(1-x)
check("[1e] Q(X) transported as a quadratic differential by X = x/(1-x) IS Q_V",
      sp.simplify(sp.together(Q_table.subs(x, g)*sp.diff(g, x)**2 - QV)) == 0)
check("NEGATIVE CONTROL omitting the Jacobian (dX/dx)^2 breaks [1e]",
      sp.simplify(sp.together(Q_table.subs(x, g) - QV)) != 0)
check("NEGATIVE CONTROL the wrong Moebius X = x/(1+x) does not give Q_V",
      sp.simplify(sp.together(Q_table.subs(x, x/(1 + x))*sp.diff(x/(1 + x), x)**2 - QV)) != 0)
check("NEGATIVE CONTROL a 1-unit change to the printed numerator breaks [1e]",
      sp.simplify(sp.together((Q_table + x**3/(4*x**2)).subs(x, g)*sp.diff(g, x)**2 - QV)) != 0)

X_sing = [sp.Integer(0), sp.Rational(1, 8), sp.Rational(1, 4),
          sp.Rational(1, 3), sp.Integer(-1), sp.oo]
pre = []
for X0 in X_sing:
    if X0 is sp.oo:
        pre.append(sp.Integer(1))                # X -> oo  <=>  x -> 1
    else:
        sol = sp.solve(sp.Eq(g, X0), x)
        pre.append(sol[0] if sol else sp.oo)     # X = -1 has no finite preimage
print("       X -> x : %s" % list(zip(X_sing, pre)))
check("[1h] singular set {0,1/8,1/4,1/3,-1,oo} corresponds to {0,1/9,1/5,1/4,oo,1}",
      pre == [sp.Integer(0), sp.Rational(1, 9), sp.Rational(1, 5),
              sp.Rational(1, 4), sp.oo, sp.Integer(1)])
print()

# ================================================================ PART 2: the Lian-Wiczer 30B entry
print("Check 2: Lian-Wiczer (arXiv:math/0611291), Conway-Norton class 30B")

QVt = QV.subs(x, sp.Symbol('t'))
t = sp.Symbol('t')
LW_num = (2700*z**8 - 2340*z**7 + 2613*z**6 + 1386*z**5
          + 311*z**4 + 112*z**3 + 15*z**2 + 2*z + 1)
LW_den_printed = 4*z**2*(60*z**4 - 37*z**3 - 25*z**2 + z + 1)**2
Q30B = LW_num/LW_den_printed

fac = sp.factor(60*z**4 - 37*z**3 - 25*z**2 + z + 1)
check("[2a] the printed denominator factors as (z-1)(3z+1)(4z+1)(5z-1)",
      sp.expand(fac - (z - 1)*(3*z + 1)*(4*z + 1)*(5*z - 1)) == 0)
print("       factor(60z^4 - 37z^3 - 25z^2 + z + 1) = %s" % fac)

t_of_z = z/(1 + 4*z)                             # 1/t = T_30B + 4, T_30B = 1/z
z_of_t = t/(1 - 4*t)
jac = sp.diff(t_of_z, z)
check("[2b] Q_V(t) dt^2 == Q_30B(z) dz^2 under t = z/(1+4z)",
      sp.simplify(sp.together(QVt.subs(t, t_of_z)*jac**2 - Q30B)) == 0)
check("[2c] the same identity read the other way, z = t/(1-4t)",
      sp.simplify(sp.together(Q30B.subs(z, z_of_t)*sp.diff(z_of_t, t)**2 - QVt)) == 0)

naive = sp.simplify(QVt.subs(t, t_of_z))
check("[2d] the naive substitution equals (1+4z)^4 * Q_30B(z): the missing Jacobian, quantified",
      sp.simplify(sp.together(naive - (1 + 4*z)**4*Q30B)) == 0)
check("NEGATIVE CONTROL the naive substitution is NOT equal to Q_30B (the spurious mismatch)",
      sp.simplify(sp.together(naive - Q30B)) != 0)
check("NEGATIVE CONTROL omitting (dz/dt)^2 in the inverse direction breaks [2c]",
      sp.simplify(sp.together(Q30B.subs(z, z_of_t) - QVt)) != 0)
check("NEGATIVE CONTROL the wrong Moebius t = z/(1+5z) fails the identity",
      sp.simplify(sp.together(QVt.subs(t, z/(1 + 5*z))*sp.diff(z/(1 + 5*z), z)**2 - Q30B)) != 0)
check("NEGATIVE CONTROL a 1-unit corruption of the tabulated numerator breaks [2b]",
      sp.simplify(sp.together(QVt.subs(t, t_of_z)*jac**2 - (LW_num + z**5)/LW_den_printed)) != 0)

poles_t = [sp.Integer(0), sp.Rational(1, 9), sp.Rational(1, 5),
           sp.Rational(1, 4), sp.Integer(1), sp.oo]
img = [sp.oo if p == sp.Rational(1, 4) else sp.nsimplify(sp.limit(z_of_t, t, p))
       for p in poles_t]
print("       t -> z : %s" % list(zip(poles_t, img)))
check("[2i] pole dictionary t={0,1/9,1/5,1/4,1,oo} -> z={0,1/5,1,oo,-1/3,-1/4}",
      img == [sp.Integer(0), sp.Rational(1, 5), sp.Integer(1), sp.oo,
              sp.Rational(-1, 3), sp.Rational(-1, 4)])

h_cusp = head(Q30B, z, 0)
h_fin = [head(Q30B, z, p) for p in
         (sp.Rational(1, 5), sp.Integer(1), sp.Rational(-1, 3), sp.Rational(-1, 4))]
Q_at_inf = sp.simplify(Q30B.subs(z, 1/s)*sp.diff(1/s, s)**2)     # local coordinate z = 1/s
h_inf = sp.simplify(sp.limit(sp.together(Q_at_inf*s**2), s, 0))
print("       heads: cusp z=0 -> %s | finite elliptic -> %s | z=oo -> %s"
      % (h_cusp, h_fin, h_inf))
check("[2j] Laurent heads 1/4 at the cusp and 3/16 at all five order-two elliptic images",
      h_cusp == sp.Rational(1, 4) and all(h == sp.Rational(3, 16) for h in h_fin)
      and h_inf == sp.Rational(3, 16))
print()

# ================================================================ PART 3: the two tables agree
print("Check 3: the two tabulated generators are related to each other")
X_of_z = sp.simplify(g.subs(x, t_of_z))          # X = x/(1-x) composed with x = t = z/(1+4z)
check("[3a] the composed Moebius map is X = z/(1+3z)",
      sp.simplify(sp.together(X_of_z - z/(1 + 3*z))) == 0)
check("[3b] Q(X) transported by that composed map IS the tabulated 30B entry",
      sp.simplify(sp.together(Q_table.subs(x, X_of_z)*sp.diff(X_of_z, z)**2 - Q30B)) == 0)
check("NEGATIVE CONTROL the composed map without its Jacobian is not the 30B entry",
      sp.simplify(sp.together(Q_table.subs(x, X_of_z) - Q30B)) != 0)
print()

print("=" * 78)
if FAILURES:
    print("RESULT: FAILURE -- %d check(s) failed:" % len(FAILURES))
    for n in FAILURES:
        print("  %s" % n)
else:
    print("RESULT: ALL CHECKS PASS")
    print("  The order-two operator underlying M is, projectively, the uniformizing equation")
    print("  tabulated for Gamma_0(30)+ in Lian-Yau and for Conway-Norton class 30B in")
    print("  Lian-Wiczer.  Both tables exhibit the entry; neither bounds it.  What the")
    print("  modularity theorem adds is the a priori pole-degree bound that makes the match")
    print("  an identity, and hence a proof of the tabulated entries.")
print("=" * 78)
sys.exit(1 if FAILURES else 0)
