"""
EXACT certificate:  L_4 = M o d/dt  is the UNIQUE factorization of L_4, and in
particular L_4 is NOT an LCLM,  L_4 != LCLM(N, d/dt)  for any order-3 operator N.

Why this matters.  An order-4 operator L with the constant solution always factors as
L = M o d/dt.  But in many lattice-model examples L turns out to be the least common left
multiple L = LCLM(N, d/dt) = N (+) d/dt instead, i.e. the constant solution is a DIRECT
SUMMAND and the order-3 operator N -- not M -- is the intrinsic object.  Everything this
paper proves about M would then have to be redone for N.  This script rules that out.

Two independent legs, each complete on its own.

  LEG A (local, at t=1).  If L_4 = LCLM(N, d/dt) then Sol(L_4) = Sol(N) (+) C, and
  d/dt : Sol(N) -> Sol(M) is an isomorphism of differential modules.  At t=1 the exponents
  of M are {-1,-1/2,1/2}; no exponent equals -1 + k for a positive integer k, so the
  Frobenius solution at the exponent -1 is a genuine single-valued Laurent series
  y = (t-1)^{-1}(1 + O(t-1)) with RESIDUE 1 != 0.  Let Y in Sol(N) be its d/dt-preimage.
  Sol(N) is monodromy-stable and meets C trivially, and (gamma Y - Y)' = gamma y - y = 0,
  so gamma Y - Y lies in Sol(N) ^ C = 0, i.e. gamma Y = Y.  But Y' = y forces
  Y = log(t-1) + (single-valued), so gamma Y = Y + 2 pi i.  Contradiction.
  CERTIFIED HERE: the exponents at t=1, the no-positive-integer-gap condition, and the
  explicit Laurent solution (exact annihilation, with a negative control).

  LEG B (global, via the adjoint).  Convention: operators compose as maps, L = A o B, and B
  is the RIGHT factor (so Sol(B) is contained in Sol(L)).  L_4 = LCLM(N, d/dt) requires an
  order-3 RIGHT factor N of L_4, i.e. L_4 = A o N with A of order one; since
  adjoint(A o N) = adjoint(N) o adjoint(A), that is equivalent to an order-1 RIGHT factor
  adjoint(A) of adjoint(L_4), equivalently a hyperexponential solution of adjoint(L_4).
  (Note L_4 = M o d/dt is the OTHER shape: there d/dt is the right factor and adjoint(M) is
  the order-3 right factor of adjoint(L_4) -- confusing the two shapes is the easy slip here.)
  Now adjoint(L_4) = (-d/dt) o adjoint(M), so such a solution y obeys adjoint(M)(y) = c, a
  constant.
    * c = 0 would make y a hyperexponential solution of adjoint(M), i.e. an order-1 right
      factor of adjoint(M) -- excluded, since M (hence adjoint(M)) is irreducible over
      Qbar(t) (certify_factor.py / CERTIFICATE.txt).
    * c != 0 forces y RATIONAL: for hyperexponential y, S := adjoint(M)(y)/y lies in
      Qbar(t), so y = c/S.
  So it suffices to show adjoint(L_4) has NO nonzero rational solution.  The space of
  rational solutions is Galois-stable, hence spanned by rational solutions over Q(t), so
  the search may be done over Q.  The search is complete: the valuation of a nonzero
  rational solution at each singular point is one of the local exponents there, which
  bounds the denominator, and the exponents at infinity bound the degree.
  CERTIFIED HERE: those bounds, and that the resulting exact linear system over Q has only
  the zero solution -- with a POSITIVE control (a manufactured L = LCLM(N, d/dt), where the
  same code does find the rational solution) and NEGATIVE controls.

Consequence (stated in the paper): the only submodules of Sol(L_4) are 0, C, and everything.
A 1-dimensional submodule maps to 0 under d/dt (its image would be a line in the irreducible
Sol(M)) hence lies in C; a 2-dimensional one would map onto a proper nonzero submodule of
Sol(M); and a 3-dimensional one is exactly the LCLM case ruled out above.  Hence
L_4 = M o d/dt is the unique factorization of L_4 into operators of positive order.

Depends on: M_coeffs.json, and (for leg B's reduction only) the irreducibility of M
certified by certify_factor.py.

Run:  python certify_lclm.py     -> writes CERTIFICATE_lclm.txt
"""
import json, os, sys
import sympy as sp

t = sp.symbols('t')
OUT = []


def say(msg):
    print(msg)
    OUT.append(msg)


FAILURES = []


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    say("  [%s] %s%s" % (tag, name, ("  " + detail) if detail else ""))
    if not ok:
        FAILURES.append(name)


# ----------------------------------------------------------------------
# Ore-algebra primitives (operator = [c0,...,cn] meaning sum ci D^i, D = d/dt)
# ----------------------------------------------------------------------
def trim(op):
    op = [sp.cancel(c) for c in op]
    while len(op) > 1 and op[-1] == 0:
        op.pop()
    return op


def Dcomp(B):
    """D o B."""
    out = [sp.Integer(0)] * (len(B) + 1)
    for k, b in enumerate(B):
        out[k] += sp.diff(b, t)
        out[k + 1] += b
    return trim(out)


def opmul(A, B):
    """A o B."""
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


def adjoint(L):
    """Formal adjoint L* = sum_k [ sum_{i>=k} (-1)^i C(i,k) a_i^{(i-k)} ] D^k."""
    n = len(L) - 1
    out = [sp.Integer(0)] * (n + 1)
    for i in range(n + 1):
        ai = L[i]
        for k in range(i + 1):
            out[k] += (-1) ** i * sp.binomial(i, k) * sp.diff(ai, t, i - k)
    return trim([sp.cancel(c) for c in out])


def apply_op(op, expr):
    """op(expr) for a sympy expression."""
    return sp.expand(sum(c * sp.diff(expr, t, i) for i, c in enumerate(op)))


def indicial_exponents(op, p):
    """Exact local exponents of `op` at the finite point t=p."""
    s, r = sp.symbols('s r')
    from collections import defaultdict
    terms = defaultdict(lambda: sp.Integer(0))
    for i, c in enumerate(op):
        ci = sp.cancel(c).subs(t, p + s)
        if ci == 0:
            continue
        coeff, m = ci.as_leading_term(s).as_coeff_exponent(s)
        fall = sp.prod([r - a for a in range(i)]) if i > 0 else sp.Integer(1)
        terms[sp.nsimplify(m) - i] += coeff * fall
    mn = min(terms.keys())
    roots = sp.roots(sp.Poly(sp.expand(terms[mn]), r))
    exps = []
    for rt, mult in roots.items():
        exps += [sp.nsimplify(rt)] * mult
    return sorted(exps, key=lambda x: sp.re(sp.N(x)))


def exponents_at_infinity(op):
    """Exponents at t=infinity, in the convention y ~ t^{-rho}."""
    s = sp.symbols('s')
    n = len(op) - 1
    expr = sum(sp.expand(c) * sp.ff(-s, i) * t ** (n - i) for i, c in enumerate(op))
    poly = sp.Poly(sp.expand(expr), t)
    ind = poly.coeff_monomial(t ** poly.degree())
    roots = sp.roots(sp.Poly(sp.expand(ind), s))
    exps = []
    for rt, mult in roots.items():
        exps += [sp.nsimplify(-rt)] * mult
    return sorted(exps, key=lambda z: sp.re(sp.N(z)))


def simple_root_exponents_mod(op, irred):
    """Exponents at a SIMPLE root of the irreducible factor `irred` of the leading
    coefficient of an order-n operator: {0,1,...,n-2} u {n-1 - c_{n-1}(r)/c_n'(r)},
    the last computed modulo `irred` and asserted to be a rational CONSTANT (so that all
    conjugate roots carry the same exponent).  Derivation: at a simple zero of c_n the
    two lowest-order Laurent contributions come from i=n and i=n-1 and the indicial
    polynomial is ff(rho,n-1) * [ c_n'(r)(rho-n+1) + c_{n-1}(r) ]."""
    n = len(op) - 1
    cn = sp.Poly(sp.expand(op[n]), t, domain='QQ')
    cn1 = sp.Poly(sp.expand(op[n - 1]), t, domain='QQ')
    ir = sp.Poly(sp.expand(irred), t, domain='QQ')
    q, rem = sp.div(cn, ir)
    assert rem.is_zero, "irred does not divide the leading coefficient"
    assert sp.gcd(ir, q).degree() == 0, "not a simple zero of the leading coefficient"
    inv = cn.diff(t).invert(ir)
    last = (sp.Poly(n - 1, t, domain='QQ') - cn1 * inv).rem(ir)
    assert last.degree() <= 0, "exponent is not constant modulo irred"
    return sorted([sp.Integer(k) for k in range(n - 1)] + [sp.nsimplify(last.as_expr())],
                  key=lambda z: sp.N(z))


def frobenius_laurent(op, p, rho, K):
    """Frobenius solution of `op` at the finite regular singular point t=p belonging to the
    exponent `rho`, as a list of K coefficients a_0=1,...,a_{K-1} of
    y = (t-p)^rho * sum_k a_k (t-p)^k.  Only valid (and only used) when no other exponent
    equals rho + k for a positive integer k, so that the recurrence never stalls."""
    s = sp.symbols('s')
    n = len(op) - 1
    csub = [sp.expand(sp.cancel(c).subs(t, p + s)) for c in op]
    # P_k(s) = s^n * sum_i c_i(p+s) * ff(rho+k, i) * s^{k-i}   (a polynomial in s)
    Pcache = {}

    def P(k):
        if k not in Pcache:
            expr = sum(csub[i] * sp.ff(rho + k, i) * s ** (n + k - i) for i in range(n + 1))
            Pcache[k] = sp.Poly(sp.expand(expr), s)
        return Pcache[k]

    def coeff(poly, m):
        return poly.coeff_monomial(s ** m) if m >= 0 else sp.Integer(0)

    a = [sp.Integer(1)]
    for m in range(1, K):
        rhs = sum(a[k] * coeff(P(k), m + n) for k in range(m))
        lead = coeff(P(m), m + n)          # = indicial polynomial at rho+m, up to scale
        assert lead != 0, "recurrence stalls at k=%d (rho+k is another exponent)" % m
        a.append(sp.cancel(-rhs / lead))
    return a


def laurent_residual(op, p, rho, a, ORD):
    """Return the list of Laurent coefficients of op(y) for the truncated
    y = (t-p)^rho sum_{k<len(a)} a_k (t-p)^k, for the orders that are FULLY determined by
    the truncation (i.e. unaffected by the discarded tail)."""
    s = sp.symbols('s')
    n = len(op) - 1
    csub = [sp.expand(sp.cancel(c).subs(t, p + s)) for c in op]
    tot = sp.Integer(0)
    for k, ak in enumerate(a):
        tot += ak * sum(csub[i] * sp.ff(rho + k, i) * s ** (n + k - i) for i in range(n + 1))
    poly = sp.Poly(sp.expand(tot), s)
    return [poly.coeff_monomial(s ** m) for m in range(0, min(ORD, len(a)) + n)]


# ----------------------------------------------------------------------
# Rational-solution search, with rigorous a-priori bounds
# ----------------------------------------------------------------------
def integer_exponents(exps):
    return [e for e in exps if sp.nsimplify(e).is_Integer]


def rational_solution_space(op, label=""):
    """Return (basis, bounds) where basis is a list of rational functions spanning the
    Q(t)-space of rational solutions of `op`.  Complete: the valuation of a nonzero
    rational solution at each singular point is a local exponent there, and the behaviour
    at infinity is governed by the exponents there."""
    n = len(op) - 1
    lead = sp.Poly(sp.expand(op[n]), t, domain='QQ')
    facs = sp.factor_list(lead.as_expr())[1]
    Q = sp.Integer(1)
    detail = []
    for f, mult in facs:
        fp = sp.Poly(sp.expand(f), t, domain='QQ')
        if fp.degree() == 0:
            continue
        if fp.degree() == 1:
            root = sp.solve(fp.as_expr(), t)[0]
            exps = indicial_exponents(op, root)
        else:
            exps = simple_root_exponents_mod(op, f)
        ints = integer_exponents(exps)
        if not ints:
            # no integer exponent => a nonzero rational solution cannot have any
            # valuation at this point at all
            return [], detail + [(str(f), exps, "no integer exponent => no rational solution")]
        d = max(0, -min(ints))
        Q *= sp.expand(f) ** d
        detail.append((str(f), exps, "pole order <= %d" % d))
    einf = exponents_at_infinity(op)
    ii = integer_exponents(einf)
    if not ii:
        return [], detail + [("infinity", einf, "no integer exponent => no rational solution")]
    # Convention-safe: whichever sign convention the indicial roots at infinity are
    # returned in, the degree of a rational solution is bounded by ONE of -min or max,
    # so taking the larger of the two is an upper bound in both cases (an over-large
    # search space can only weaken a POSITIVE find, never a negative conclusion).
    dinf = max(-min(ii), max(ii))
    detail.append(("infinity", einf, "deg(y) <= %d" % dinf))
    Q = sp.expand(Q)
    degQ = sp.Poly(Q, t).degree() if Q != 1 else 0
    degP = degQ + int(dinf)
    if degP < 0:
        return [], detail
    # Gauge the operator by 1/Q so the unknown is a POLYNOMIAL, and clear denominators:
    #   op(z/Q) = 0   <=>   sum_k b_k(t) z^{(k)} = 0,
    #   b_k = Q^{n+1} * sum_{i>=k} C(i,k) a_i (1/Q)^{(i-k)}      (all b_k polynomial).
    invd = []                                   # invd[j] = Q^{n+1} * (1/Q)^{(j)}
    cur = sp.Integer(1) / Q
    for j in range(n + 1):
        invd.append(sp.expand(sp.cancel(cur * Q ** (n + 1))))
        cur = sp.diff(cur, t)
    b = []
    for k in range(n + 1):
        bk = sum(sp.binomial(i, k) * sp.expand(op[i]) * invd[i - k] for i in range(k, n + 1))
        b.append(sp.expand(bk))
    # columns of the linear system: the image of each monomial t^j
    cols = []
    for j in range(degP + 1):
        img = sp.expand(sum(b[k] * sp.diff(t ** j, t, k) for k in range(n + 1)))
        cols.append(sp.Poly(img, t) if img != 0 else sp.Poly(0, t))
    nzcols = [c for c in cols if c.as_expr() != 0]
    if not nzcols:
        Mat = sp.zeros(1, degP + 1)
    else:
        maxdeg = max(c.degree() for c in nzcols)
        Mat = sp.Matrix([[c.coeff_monomial(t ** d) for c in cols]
                         for d in range(maxdeg + 1)])
    ns = Mat.nullspace()
    basis = []
    for vec in ns:
        Pv = sp.expand(sum(vec[j] * t ** j for j in range(degP + 1)))
        if Pv != 0:
            basis.append(sp.cancel(Pv / Q))
    return basis, detail


def operator_from_solutions(sols):
    """Monic operator (as a list) annihilating exactly the given independent functions,
    via the Wronskian determinant expansion."""
    n = len(sols)
    rows = [[sp.diff(f, t, i) for f in sols] for i in range(n + 1)]
    coeffs = []
    for i in range(n + 1):
        sub = sp.Matrix([rows[r] for r in range(n + 1) if r != i])
        coeffs.append((-1) ** (n - i) * sp.cancel(sub.det()))
    lead = coeffs[n]
    return trim([sp.cancel(c / lead) for c in coeffs])


# ======================================================================
say("=" * 74)
say("certify_lclm.py -- L_4 = M o d/dt is the UNIQUE factorization; no LCLM")
say("=" * 74)

# ----------------------------------------------------------------------
# 0. VALIDATION of the primitives on operators of KNOWN structure
# ----------------------------------------------------------------------
say("")
say("[0] validating primitives on operators of known structure")

# adjoint is an involution
Ltest = [t ** 2 + 1, t, sp.Integer(3), sp.Integer(1)]
check("adjoint involutive", all(sp.cancel(a - b) == 0
                                for a, b in zip(Ltest, adjoint(adjoint(Ltest)))))

# adjoint(A o B) = adjoint(B) o adjoint(A)
Aop = [sp.Integer(1), t, sp.Integer(1)]
Bop = [-1 / t, sp.Integer(1)]
lhs = adjoint(opmul(Aop, Bop))
rhs = opmul(adjoint(Bop), adjoint(Aop))
check("adjoint(A o B) = adjoint(B) o adjoint(A)",
      len(lhs) == len(rhs) and all(sp.cancel(x - y) == 0 for x, y in zip(lhs, rhs)))

# operator_from_solutions reproduces a known operator
Dfour = operator_from_solutions([sp.Integer(1), t, t ** 2, t ** 3])
check("operator_from_solutions({1,t,t^2,t^3}) = D^4",
      trim(Dfour) == [sp.Integer(0)] * 4 + [sp.Integer(1)])

# Frobenius machinery on Euler operators with known exponents.
# (a) exponents {-1,1/2}: NO positive-integer gap above -1, exact solutions t^-1, t^{1/2}.
EulerOK = [sp.Rational(-1, 2), sp.Rational(3, 2) * t, t ** 2]
check("Euler control (a) exponents {-1,1/2}",
      indicial_exponents(EulerOK, 0) == [sp.Rational(-1), sp.Rational(1, 2)])
aE = frobenius_laurent(EulerOK, 0, sp.Integer(-1), 6)
check("Euler control (a): Frobenius at -1 returns t^{-1} exactly",
      aE[0] == 1 and all(x == 0 for x in aE[1:]))
# (b) exponents {-1,2}: -1 + 3 IS another exponent, so the log-free construction must
#     REFUSE to run.  This is the guard that makes the t=1 gap test non-vacuous.
EulerStall = [sp.Integer(-2), sp.Integer(0), t ** 2]
check("Euler control (b) exponents {-1,2}", indicial_exponents(EulerStall, 0) == [-1, 2])
stalled = False
try:
    frobenius_laurent(EulerStall, 0, sp.Integer(-1), 6)
except AssertionError:
    stalled = True
check("Euler control (b): the log-free construction correctly refuses to run", stalled)

# rational-solution solver: positive control with a KNOWN pole
L1 = [sp.Integer(2) / (t - 2), sp.Integer(1)]           # y' + 2y/(t-2) = 0  -> y = (t-2)^{-2}
L1 = [sp.expand(c * (t - 2)) for c in L1]               # clear denominators
b1, _ = rational_solution_space(L1)
check("solver finds the planted pole solution (t-2)^{-2}",
      len(b1) == 1 and sp.cancel(b1[0] * (t - 2) ** 2).is_number)

# rational-solution solver: negative control (no rational solution)
L2 = [sp.Integer(-1), t ** 2]                            # t^2 y' = y -> y = exp(-1/t)
b2, _ = rational_solution_space(L2)
check("solver returns nothing for exp(-1/t) (negative control)", b2 == [])

say("")
say("[0b] end-to-end POSITIVE control: a manufactured L = LCLM(N, d/dt)")
Lsplit = operator_from_solutions([sp.Integer(1), t, t ** 2, sp.Integer(1) / (t - 2)])
den = sp.lcm([sp.denom(sp.cancel(c)) for c in Lsplit])
Lsplit = trim([sp.expand(sp.cancel(c * den)) for c in Lsplit])
check("control L has the constant solution", sp.cancel(apply_op(Lsplit, sp.Integer(1))) == 0)
bsplit, _ = rational_solution_space(adjoint(Lsplit))
check("adjoint(control) HAS a nonzero rational solution => split detected",
      len(bsplit) >= 1)

# ----------------------------------------------------------------------
# 1. Load M and build L_4
# ----------------------------------------------------------------------
say("")
say("[1] the hyperkagome operator")
HERE = os.path.dirname(os.path.abspath(__file__))
Mc = json.load(open(os.path.join(HERE, "M_coeffs.json")))
assert Mc["order"] == 3
M = [sum(int(c) * t ** j for j, c in enumerate(Mc["coeffs"][i])) for i in range(4)]
L4 = opmul(M, [sp.Integer(0), sp.Integer(1)])            # L_4 = M o d/dt
say("  M order 3, deg c3 = %d" % sp.Poly(M[3], t).degree())
check("L_4 = M o d/dt has order 4", len(L4) - 1 == 4)
check("L_4 annihilates the constant", sp.cancel(apply_op(L4, sp.Integer(1))) == 0)
check("L_4 leading coefficient equals that of M", sp.expand(L4[4] - M[3]) == 0)

# ----------------------------------------------------------------------
# 2. LEG A -- the residue obstruction at t=1
# ----------------------------------------------------------------------
say("")
say("[2] LEG A: the exponent -1 solution of M at t=1 has a nonzero residue")
e1 = indicial_exponents(M, sp.Integer(1))
say("  exponents of M at t=1: %s" % [str(x) for x in e1])
check("exponents at t=1 are {-1,-1/2,1/2}",
      sorted(e1, key=lambda z: sp.N(z)) == [sp.Rational(-1), sp.Rational(-1, 2), sp.Rational(1, 2)])
check("the exponent -1 is a SIMPLE root of the indicial equation at t=1",
      list(e1).count(sp.Rational(-1)) == 1 and len(e1) == 3)
gaps = [e - sp.Integer(-1) for e in e1]
bad = [g for g in gaps if sp.nsimplify(g).is_Integer and g > 0]
check("no exponent equals -1 + (positive integer) => the -1 solution is log-free",
      bad == [], "gaps %s" % [str(g) for g in gaps])

KTERMS = 30
aL = frobenius_laurent(M, sp.Integer(1), sp.Integer(-1), KTERMS)
check("leading Laurent coefficient (= the residue at t=1) is 1", aL[0] == 1)
resid = laurent_residual(M, sp.Integer(1), sp.Integer(-1), aL, KTERMS)
nz = [m for m, c in enumerate(resid) if sp.cancel(c) != 0]
check("M annihilates the Laurent solution exactly through order %d" % (KTERMS - 1),
      all(sp.cancel(c) == 0 for c in resid[:KTERMS]),
      "first nonzero residual order: %s" % (nz[0] if nz else "none in range"))

# NEGATIVE CONTROL: corrupting one coefficient must break the annihilation
aBad = list(aL)
aBad[3] = aBad[3] + 1
residBad = laurent_residual(M, sp.Integer(1), sp.Integer(-1), aBad, KTERMS)
check("negative control: corrupting a_3 breaks the annihilation",
      any(sp.cancel(c) != 0 for c in residBad[:KTERMS]))

say("  => any antiderivative of that solution is log(t-1) + (single-valued),")
say("     so no monodromy-stable complement of the constants can contain it.")

# ----------------------------------------------------------------------
# 3. LEG B -- adjoint(L_4) has no rational solution
# ----------------------------------------------------------------------
say("")
say("[3] LEG B: adjoint(L_4) has NO nonzero rational solution")
A4 = adjoint(L4)
check("adjoint(L_4) has order 4", len(A4) - 1 == 4)
# the structural identity the reduction rests on: adjoint(M o d/dt) = (-d/dt) o adjoint(M),
# so a solution y of adjoint(L_4) is exactly one with adjoint(M)(y) = constant.
AM = adjoint(M)
comp = opmul([sp.Integer(0), sp.Integer(-1)], AM)
check("adjoint(L_4) = (-d/dt) o adjoint(M)",
      len(comp) == len(A4) and all(sp.expand(x - y) == 0 for x, y in zip(comp, A4)))
# cross-check the mod-irreducible exponent formula (used on p_7) against the direct
# indicial computation at a genuine SIMPLE rational root of the same leading coefficient.
e_direct = indicial_exponents(A4, sp.Rational(1, 4))
e_mod = simple_root_exponents_mod(A4, 4 * t - 1)
check("simple-root exponent formula agrees with the direct indicial computation at t=1/4",
      sorted(e_direct, key=lambda z: sp.N(z)) == sorted(e_mod, key=lambda z: sp.N(z)),
      "direct %s vs mod %s" % ([str(x) for x in e_direct], [str(x) for x in e_mod]))
basis, detail = rational_solution_space(A4)
for f, exps, note in detail:
    say("  at %-28s exponents %-34s %s"
        % (f, "{" + ",".join(str(x) for x in exps) + "}", note))
check("no nonzero rational solution of adjoint(L_4)", basis == [],
      "" if basis == [] else "found %s" % basis)

# NEGATIVE CONTROL: the same solver on M itself (irreducible => no rational solution)
bM, _ = rational_solution_space(M)
check("negative control: solver returns nothing for the irreducible M", bM == [])

say("")
say("  Chain: a rational solution of adjoint(L_4) is the only possible hyperexponential")
say("  one (adjoint(M)(y) = c; c = 0 would contradict the certified irreducibility of M,")
say("  and c != 0 forces y = c/S rational).  None exists, so adjoint(L_4) has no order-1")
say("  right factor, so L_4 has no order-3 right factor.")

# ----------------------------------------------------------------------
# 4. CERTIFICATE
# ----------------------------------------------------------------------
say("")
say("=" * 74)
if FAILURES:
    say("CERTIFICATE INCOMPLETE -- failures: %s" % ", ".join(FAILURES))
else:
    say("RESULT: L_4 has no order-3 right factor, hence L_4 != LCLM(N, d/dt) for every")
    say("order-3 operator N: the constant solution is NOT a direct summand.  With M")
    say("irreducible (certify_factor.py), the only submodules of Sol(L_4) are 0, the")
    say("constants, and everything, so L_4 = M o d/dt is the UNIQUE factorization of L_4")
    say("into operators of positive order.")
    say("Proved twice, independently: LEG A (residue of the exponent -1 solution at t=1)")
    say("and LEG B (no rational solution of adjoint(L_4), with complete a-priori bounds).")
say("=" * 74)

with open(os.path.join(HERE, "CERTIFICATE_lclm.txt"), "w") as f:
    f.write("\n".join(OUT) + "\n")
print("\nwrote CERTIFICATE_lclm.txt")
sys.exit(1 if FAILURES else 0)
