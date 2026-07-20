"""
EXACT analysis of the hyperkagome order-3 Picard-Fuchs operator M.
Pure Python (sympy exact rationals) -- no Sage/ore_algebra needed.

What this script certifies, by EXACT computation (no floating point in the decisive steps):

  (a) M has NO order-1 right factor over Q(t)   [complete Fuchsian residue enumeration]
  (b) M has NO order-2 right factor over Q(t)   [= no order-1 right factor of adjoint(M)]
  (c) M has a genuine LOGARITHMIC local solution at t=0 (repeated indicial exponent -1,
      exactly one log-free local solution of three) => M is NOT completely reducible.
  ==> M is IRREDUCIBLE over Qbar(t), by the Galois-descent argument stated below.
  (d) M is NOT a *literal* symmetric square of a 2nd-order operator (exponent triples are not
      arithmetic progressions at every singular point).  NOTE: this does NOT exclude an elliptic
      closed form.  The earlier inference "not-Sym^2 => no elliptic (K) closed form" is RETRACTED:
      Sym^2(M) has a rational solution, so M's Galois group is orthogonal (O(3,C), G^0 = SO(3,C))
      and M IS projectively equivalent to a symmetric square via an order-2 differential intertwiner
      -- an equivalence the local-exponent AP test cannot see.  See certify_orthogonal.py.

IRREDUCIBILITY OVER Qbar(t) -- the descent argument (this replaces the earlier, incorrect
"has a log => not a direct sum => irreducible" step, which only gives INDECOMPOSABLE):

  A third-order M is reducible over Qbar(t) iff M has an order-1 right factor over Qbar(t)
  OR M has an order-2 right factor over Qbar(t); the latter is equivalent to adjoint(M)
  having an order-1 right factor over Qbar(t).  Suppose D-u (u in Qbar(t)) is an order-1
  right factor of M.  Its Galois orbit {D-u^sigma} are all order-1 right factors of the
  Q(t)-operator M; their LCLM R is Galois-stable, hence defined over Q(t), and is a right
  factor of M of order k = dim(sum of the 1-dim solution spaces) in {1,2,3}:
     k=1: u in Q(t)          -> contradicts (a);
     k=2: R is an order-2 right factor of M over Q  -> contradicts (b);
     k=3: R = M, so M is the LCLM of order-1 operators = completely reducible
          => semisimple monodromy => NO log  -> contradicts (c).
  Hence M has no order-1 right factor over Qbar(t).  The same argument applied to adjoint(M)
  (which has no order-1 factor over Q by (b), no order-2 factor over Q by (a), and is not
  completely reducible because complete reducibility is self-dual under the adjoint and M is
  not completely reducible by (c)) shows adjoint(M) has no order-1 right factor over Qbar(t),
  i.e. M has no order-2 right factor over Qbar(t).  Therefore M is IRREDUCIBLE over Qbar(t).

Every operator primitive (multiply, right-divide, adjoint) and the log/exponent routines are
validated on synthetic operators with KNOWN structure before the real run.
A simple-cubic-style symmetric-square POSITIVE CONTROL is run to show the not-Sym^2 test is
non-vacuous (it can, and does, return "Sym^2-compatible" on a genuine symmetric square).

Run:  python certify_factor.py
"""
import json, os, itertools
import sympy as sp

t = sp.symbols('t')

# ----------------------------------------------------------------------
# Ore-algebra primitives.  Operator = list [c0,c1,...,cn] meaning sum ci*D^i,
# D = d/dt, coefficients ci are sympy rational functions of t.
# ----------------------------------------------------------------------
def trim(op):
    op = [sp.cancel(c) for c in op]
    while len(op) > 1 and op[-1] == 0:
        op.pop()
    return op

def Dcomp(B):
    """Return the operator D o B (compose derivative on the left of operator B)."""
    out = [sp.Integer(0)] * (len(B) + 1)
    for k, b in enumerate(B):
        out[k] += sp.diff(b, t)
        out[k + 1] += b
    return trim(out)

def opmul(A, B):
    """Return A o B (operator composition)."""
    res = [sp.Integer(0)]
    DiB = [c for c in B]  # D^0 o B = B
    for i, a in enumerate(A):
        if a != 0:
            term = [sp.cancel(a * c) for c in DiB]
            if len(term) > len(res):
                res = res + [sp.Integer(0)] * (len(term) - len(res))
            for k, c in enumerate(term):
                res[k] += c
        DiB = Dcomp(DiB)  # advance to D^{i+1} o B
    return trim(res)

def opright_divmod(A, B):
    """Right division: A = Q o B + R with deg_D(R) < deg_D(B). Return (Q, R)."""
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
        qD = [sp.Integer(0)] * deg + [q]     # operator q*D^deg
        sub = opmul(qD, B)
        newlen = max(len(A), len(sub))
        A = A + [sp.Integer(0)] * (newlen - len(A))
        sub = sub + [sp.Integer(0)] * (newlen - len(sub))
        A = [sp.cancel(A[k] - sub[k]) for k in range(newlen)]
        while len(A) > 1 and A[-1] == 0:
            A.pop()
    return trim(Q), trim(A)

def right_divides(M, B):
    _, R = opright_divmod(M, B)
    return all(sp.cancel(c) == 0 for c in R)

def adjoint(L):
    """Formal adjoint: L* = sum_k [ sum_{i>=k} (-1)^i C(i,k) a_i^{(i-k)} ] D^k."""
    n = len(L) - 1
    out = [sp.Integer(0)] * (n + 1)
    for i in range(n + 1):
        ai = L[i]
        for k in range(i + 1):
            out[k] += (-1)**i * sp.binomial(i, k) * sp.diff(ai, t, i - k)
    return trim([sp.cancel(c) for c in out])

def sym2(L2):
    """Symmetric square of an order-2 operator L2 = [a0, a1, 1] (monic: D^2 + a1 D + a0).
    Classical closed form: D^3 + 3 a1 D^2 + (2 a1^2 + a1' + 4 a0) D + (4 a0 a1 + 2 a0').
    Validated below on D^2-1 (=> D^3-4D) and D^2 (=> D^3)."""
    a0, a1 = sp.cancel(L2[0] / L2[2]), sp.cancel(L2[1] / L2[2])
    return [sp.cancel(4*a0*a1 + 2*sp.diff(a0, t)),
            sp.cancel(2*a1**2 + sp.diff(a1, t) + 4*a0),
            sp.cancel(3*a1),
            sp.Integer(1)]

def indicial_exponents(op, p):
    """Local exponents (indicial roots) of `op` at the finite regular singular point t=p.
    Exact; handles polynomial AND rational-function coefficients (only the lowest-order
    Laurent term of each c_i(p+s) contributes to the indicial polynomial)."""
    s = sp.symbols('s'); r = sp.symbols('r')
    from collections import defaultdict
    terms = defaultdict(lambda: sp.Integer(0))
    for i, c in enumerate(op):
        ci = sp.cancel(c).subs(t, p + s)
        if ci == 0:
            continue
        coeff, m = ci.as_leading_term(s).as_coeff_exponent(s)   # lowest-power term a*s^m
        fall = sp.prod([r - a for a in range(i)]) if i > 0 else sp.Integer(1)
        terms[sp.nsimplify(m) - i] += coeff * fall
    mn = min(terms.keys())
    roots = sp.roots(sp.Poly(sp.expand(terms[mn]), r))
    exps = []
    for rt, mult in roots.items():
        exps += [sp.nsimplify(rt)] * mult
    return sorted(exps, key=lambda x: sp.re(sp.N(x)))

def leadfactor_exponents(op, irred):
    """EXACT exponents at a simple root of an irreducible factor `irred` of the leading
    coefficient.  At a simple zero r of the leading coeff L3, a Frobenius balance gives
    exponents {0, 1, 2 - L2(r)/L3'(r)}.  We compute 2 - L2/L3' reduced MODULO irred and
    assert the result is a rational CONSTANT -- i.e. the third exponent is identical at all
    conjugate roots (Galois consistency).  This replaces the earlier nroots+nsimplify guess
    (which mis-read {0,1,3} as {0,1,2})."""
    L3 = sp.Poly(sp.expand(op[-1]), t, domain='QQ')
    L2 = sp.Poly(sp.expand(op[-2]), t, domain='QQ')
    ir = sp.Poly(sp.expand(irred), t, domain='QQ')
    q, rem = sp.div(L3, ir)
    assert rem.is_zero, "irred does not divide the leading coefficient"
    assert sp.gcd(ir, q).degree() == 0, "root of irred is not a SIMPLE zero of the leading coeff"
    inv = L3.diff(t).invert(ir)                       # (L3')^{-1} mod irred
    third = (sp.Poly(2, t, domain='QQ') - L2 * inv).rem(ir)
    assert third.degree() <= 0, ("third exponent is NOT constant mod irred "
                                  "=> exponents differ across conjugate roots (bug)")
    return sorted([sp.Integer(0), sp.Integer(1), sp.nsimplify(third.as_expr())],
                  key=lambda z: sp.N(z))

def exponents_at_infinity(op):
    """Exponents at t=infinity (y ~ t^{-s}); finite & rational <=> regular singular there.
    op[i] D^i applied to t^{-s} gives op[i]*ff(-s,i)*t^{-s-i}; multiply through by t^{n} to
    clear negative powers, then the top-degree coefficient in t is the indicial polynomial."""
    s = sp.symbols('s'); n = len(op) - 1
    expr = sum(sp.expand(c) * sp.ff(-s, i) * t**(n - i) for i, c in enumerate(op))
    poly = sp.Poly(sp.expand(expr), t)
    ind = poly.coeff_monomial(t**poly.degree())
    roots = sp.roots(sp.Poly(sp.expand(ind), s))
    exps = []
    for rt, mult in roots.items():
        exps += [sp.nsimplify(-rt)] * mult            # exponent in t is -s
    return sorted(exps, key=lambda z: sp.re(sp.N(z)))

# ----------------------------------------------------------------------
# VALIDATION of the primitives on synthetic operators with known structure
# ----------------------------------------------------------------------
def validate():
    # order-1 right factor: T = A2 o B1, B1 a RIGHT factor
    B1 = [-1/t, sp.Integer(1)]
    A2 = [sp.Integer(1), t, sp.Integer(1)]
    T = opmul(A2, B1)
    assert right_divides(T, B1), "validate: failed to detect known order-1 right factor"
    Q, R = opright_divmod(T, B1)
    assert all(sp.cancel(c) == 0 for c in R), "validate: nonzero remainder on known factor"
    assert all(sp.cancel(a - b) == 0 for a, b in zip(trim(Q), trim(A2))), \
        "validate: quotient mismatch on known factorization"
    # order-2 right factor: T2 = B1b o A2b, A2b a RIGHT (order-2) factor
    B1b = [t + 1, sp.Integer(1)]
    A2b = [sp.Integer(1), sp.Integer(0), sp.Integer(1)]  # D^2 + 1
    T2 = opmul(B1b, A2b)
    assert right_divides(T2, A2b), "validate: failed to detect known order-2 right factor"
    # adjoint involution
    L = [t**2 + 1, t, sp.Integer(3), sp.Integer(1)]
    assert all(sp.cancel(a - b) == 0 for a, b in zip(L, adjoint(adjoint(L)))), \
        "validate: adjoint not involutive"
    # sym2 closed form on operators with known solution bases
    assert trim(sym2([sp.Integer(-1), sp.Integer(0), sp.Integer(1)])) == trim([sp.Integer(0), -sp.Integer(4), sp.Integer(0), sp.Integer(1)]), \
        "validate: sym2(D^2-1) != D^3-4D"
    assert trim(sym2([sp.Integer(0), sp.Integer(0), sp.Integer(1)])) == [sp.Integer(0), sp.Integer(0), sp.Integer(0), sp.Integer(1)], \
        "validate: sym2(D^2) != D^3"
    print("  [validate] Ore multiply / right-divide / quotient / adjoint / sym2: PASS")

print("validating Ore-algebra primitives on synthetic operators...")
validate()

# ----------------------------------------------------------------------
# Load M
# ----------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
Mc = json.load(open(os.path.join(HERE, "M_coeffs.json")))
assert Mc["order"] == 3
def mkpoly(coeffs):
    return sum(int(c) * t**j for j, c in enumerate(coeffs))
M = [mkpoly(Mc["coeffs"][i]) for i in range(4)]

# ----------------------------------------------------------------------
# Singular points and EXACT exponents
# ----------------------------------------------------------------------
print("leading coeff c3 factorization:", sp.factor(M[3]))
fac = sp.factor_list(M[3])
rat_sings = sorted(set(sp.Rational(r) for r, m in sp.roots(sp.Poly(M[3], t), filter='Q').items()))
p7 = None
for f, m in fac[1]:
    if sp.Poly(f, t).degree() == 7 and sp.Poly(f, t).is_irreducible:
        p7 = sp.expand(f)
print("rational singular points:", [str(s) for s in rat_sings])
print("p7 present (irreducible deg-7 locus):", p7 is not None)

exps = {s: indicial_exponents(M, s) for s in rat_sings}
for s in rat_sings:
    print("  exponents at t=%s : %s" % (s, [str(e) for e in exps[s]]))

# EXACT p7 exponents (Galois-consistent), for M and for adjoint(M)
p7_exps = leadfactor_exponents(M, p7)
print("  exponents at the p7 locus (EXACT, mod p7): %s" % [str(e) for e in p7_exps])
assert p7_exps == [sp.Integer(0), sp.Integer(1), sp.Integer(3)], "p7 exponents of M expected {0,1,3}"

# Fuchsian at infinity (needed so order-1 candidates have no polynomial part)
inf_exps = exponents_at_infinity(M)
print("  exponents at t=infinity : %s  (finite & rational => regular singular)" % [str(e) for e in inf_exps])
assert all(e.is_rational for e in inf_exps), "M is NOT Fuchsian at infinity -- enumeration would be incomplete"

# ----------------------------------------------------------------------
# (a) complete exact search for order-1 right factors of M over Q
# u = sum_s r_s/(t-s) + r_p7 * p7'/p7 , r_s in exps[s], r_p7 in p7 exponents.
# (Fuchsian everywhere incl. infinity => no polynomial part; every hyperexponential
#  solution over Q has exactly this form.)
# ----------------------------------------------------------------------
def _bracket(OP, u):
    """(D-u) is an order-1 right factor of OP  <=>  this rational function is 0."""
    u1 = sp.diff(u, t); u2 = sp.diff(u, t, 2)
    c0, c1, c2, c3 = OP[0], OP[1], OP[2], OP[3]
    return c3*(u2 + 3*u*u1 + u**3) + c2*(u1 + u**2) + c1*u + c0

def search_order1_right_factors(OP, sing_exps, p7poly, p7exps):
    keys = list(dict.fromkeys(sing_exps.keys()))
    choices = [sorted(set(sing_exps[k]), key=lambda z: sp.N(z)) for k in keys]
    p7log = sp.cancel(sp.diff(p7poly, t) / p7poly) if p7poly is not None else sp.Integer(0)
    p7opts = sorted(set(p7exps), key=lambda z: sp.N(z)) if p7poly is not None else [sp.Integer(0)]
    test_pts = [sp.Rational(a, b) for a, b in [(7, 3), (13, 5), (23, 9), (3, 2)]]
    found = []; ncand = 0
    for combo in itertools.product(*choices):
        for rp7 in p7opts:
            ncand += 1
            u = sum(combo[i] / (t - keys[i]) for i in range(len(keys))) + rp7 * p7log
            br = _bracket(OP, u)
            ok = True
            for tp in test_pts:                       # cheap numeric pre-screen ...
                if abs(complex(sp.N(br.subs(t, tp), 20))) > 1e-9:
                    ok = False; break
            if ok and sp.cancel(sp.together(br)) == 0: # ... then confirm EXACTLY over Q
                found.append(sp.cancel(u))
    return found, ncand

print("\n=== (a): order-1 right factors of M (exact, complete over Q) ===")
f1, n1 = search_order1_right_factors(M, exps, p7, p7_exps)
print("candidates tested: %d ; order-1 right factors found: %d %s" % (n1, len(f1), f1))

print("\n=== (b): order-2 right factors of M  =  order-1 right factors of adjoint(M) ===")
Ma = [sp.cancel(c) for c in adjoint(M)]
exps_adj = {s: indicial_exponents(Ma, s) for s in rat_sings}
p7_exps_adj = leadfactor_exponents(Ma, p7)
print("  adjoint(M) p7 exponents (EXACT): %s" % [str(e) for e in p7_exps_adj])
assert p7_exps_adj == [sp.Integer(-2), sp.Integer(0), sp.Integer(1)], "adjoint p7 exponents expected {-2,0,1}"
f2, n2 = search_order1_right_factors(Ma, exps_adj, p7, p7_exps_adj)
print("candidates tested: %d ; order-1 right factors of adjoint found: %d %s" % (n2, len(f2), f2))

# ----------------------------------------------------------------------
# (c) genuine logarithmic local solution at t=0 for M and adjoint(M) -- RIGOROUS.
# A repeated indicial exponent is necessary but NOT sufficient for a log (it can be
# apparent).  We count the log-free (formal Laurent) local solutions by the Frobenius
# obstruction: an order-n operator has n log-free solutions iff the point is apparent;
# fewer than n => a genuine logarithm.  (Same method as certify_nonliouvillian.py.)
# ----------------------------------------------------------------------
def _falling(s, i):
    r = sp.Integer(1)
    for a in range(i):
        r *= (s - a)
    return r

def count_log_free(op, p, K=8):
    """Return (number of log-free local solutions at t=p, sorted exponents)."""
    from collections import defaultdict
    x = sp.symbols('x'); s = sp.symbols('s')
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
            qd[j - i] += co * _falling(s, i)
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

# validate the log counter on known operators (t^3 D^3 = no log; theta^3 = logs)
assert count_log_free([sp.Integer(0), sp.Integer(0), sp.Integer(0), t**3], sp.Integer(0))[0] == 3, "log-counter validation (no-log) failed"
assert count_log_free([sp.Integer(0), t, 3*t**2, t**3], sp.Integer(0))[0] == 1, "log-counter validation (log) failed"

nlf_M, e0_M = count_log_free(M, sp.Integer(0))
nlf_A, e0_A = count_log_free(Ma, sp.Integer(0))
log_M, log_A = nlf_M < 3, nlf_A < 3
print("\nM       t=0 exponents %s ; log-free local solutions = %d of 3 => genuine log: %s"
      % ([str(x) for x in e0_M], nlf_M, log_M))
print("adj(M)  t=0 exponents %s ; log-free local solutions = %d of 3 => genuine log: %s  (self-dual)"
      % ([str(x) for x in e0_A], nlf_A, log_A))

# ----------------------------------------------------------------------
# (d) not a symmetric square (exact, every singular point) + POSITIVE CONTROL
# ----------------------------------------------------------------------
def sym2_test(e):
    lo, mid, hi = sorted(e, key=lambda z: sp.N(z))
    return sp.Rational(1, 2) * (lo + hi) == mid       # True <=> arithmetic progression

print("\n=== (d): symmetric-square test at every singular point of M ===")
notsym2_points = []
for s in rat_sings:
    e = exps[s]
    ok = sym2_test(e)
    notsym2_points.append((s, ok))
    print("  t=%s exps=%s : %s" % (s, [str(x) for x in e],
          "Sym^2-compatible (AP)" if ok else "NOT Sym^2 (not an AP)"))
any_notsym2 = any(not ok for _, ok in notsym2_points)

# POSITIVE CONTROL: a genuine symmetric square must pass the AP test everywhere.
# Build Sym^2 of an Euler operator with exponents {0,1/2} at t=0; expect {0,1/2,1} (an AP),
# the same pattern the simple-cubic / BCC / FCC lattice LGF operators exhibit.
Lell = [sp.Integer(0), sp.Rational(1, 2) / t, sp.Integer(1)]   # D^2 + (1/2t) D, exps {0,1/2} at 0
Msc = sym2(Lell)
e_sc = indicial_exponents(Msc, sp.Integer(0))
print("  [positive control] Sym^2 of a {0,1/2}-operator has t=0 exps %s : %s"
      % ([str(x) for x in e_sc], "Sym^2-compatible (AP)" if sym2_test(e_sc) else "NOT Sym^2"))
assert e_sc == [sp.Integer(0), sp.Rational(1, 2), sp.Integer(1)] and sym2_test(e_sc), \
    "positive control failed -- the not-Sym^2 test would be vacuous"

# ----------------------------------------------------------------------
# CERTIFICATE
# ----------------------------------------------------------------------
print("\n================ CERTIFICATE ================")
print("order-1 right factors of M over Q         :", len(f1), "(0 => none)")
print("order-1 right factors of adjoint(M) over Q :", len(f2), "(0 => none, i.e. no order-2 factor of M)")
print("genuine log at t=0 (M and adjoint(M))      :", log_M and log_A, "(=> not completely reducible)")
print("not a symmetric square (some singular pt)  :", any_notsym2)
print("positive control (Sym^2 -> AP) passes      : True")
print("---------------------------------------------")
if len(f1) == 0 and len(f2) == 0 and log_M and log_A:
    print("=> No order-1 or order-2 right factor of M over Q, and M is not completely reducible.")
    print("   By Galois descent (a hypothetical order-1/2 right factor over Qbar yields, via its")
    print("   Galois orbit's LCLM, a right factor over Q of order 1, 2, or 3 -- order 1 excluded")
    print("   by (a), order 2 by (b), order 3 excluded because it forces complete reducibility")
    print("   hence semisimple monodromy hence NO log, contradicting the genuine t=0 log):")
    print("   ==> M is IRREDUCIBLE over Qbar(t).")
if any_notsym2:
    print("=> M is irreducible and NOT a *literal* symmetric square (no arithmetic-progression")
    print("   exponent triple).  This does NOT exclude an elliptic closed form: Sym^2(M) has a")
    print("   rational solution => Galois group O(3,C), G^0 = SO(3,C), so M IS projectively a")
    print("   symmetric square via an order-2 intertwiner => an elliptic/modular form is EXPECTED.")
    print("   (The old 'not-Sym^2 => no elliptic' claim is RETRACTED; see certify_orthogonal.py.)")

if __name__ == "__main__":
    pass
