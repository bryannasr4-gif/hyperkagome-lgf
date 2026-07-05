"""
STRENGTHENING #2: show that the hyperkagome operator M has NO Liouvillian solutions, i.e. no
ALGEBRAIC or ELEMENTARY closed form.  This is a SEPARATE exclusion from "no elliptic closed
form" (which comes from not-being-a-symmetric-square, in certify_factor.py).  Combining the
two legs gives: the hyperkagome LGF has no closed form in algebraic, elementary, OR elliptic
functions.  NOTE: non-Liouvillian does NOT by itself exclude an elliptic form -- the complete
elliptic integral K is itself non-Liouvillian, and the simple-cubic LGF operator is
irreducible AND non-Liouvillian yet has a closed form in K.  So the elliptic exclusion must
come from the independent not-Sym^2 argument, and only the two together give the full result.
This whole conclusion is conditional on M being the minimal annihilator of the LGF (the
guess-and-verify standard; an unconditional creative-telescoping proof remains open).

Argument (all pieces exact):
  Let G be the differential Galois group of the irreducible order-3 operator M over Qbar(t).
  By Singer-Ulmer, an irreducible order-3 operator has Liouvillian solutions  <=>  G is
  FINITE or IMPRIMITIVE.
    * FINITE  => every solution algebraic => no logarithms anywhere.
    * IMPRIMITIVE => G permutes 3 lines, so in the block basis every group element (in
      particular every local monodromy) is a MONOMIAL matrix; an invertible monomial matrix
      is semisimple (some power is diagonal), so its Jordan form is diagonal => NO logarithm.
  Hence: if M has a genuine LOGARITHMIC local solution, G is neither finite nor imprimitive,
  so G acts irreducibly and non-solvably => M is NON-LIOUVILLIAN (no algebraic/elementary
  solutions).  [The not-Sym^2 fact separately excludes the elliptic/K case.]

This script PROVES a logarithm exists at t=0 by counting log-free (formal Laurent-series)
local solutions: an order-3 operator has 3 log-free solutions iff the point is "apparent"
(no log). We show M has < 3 log-free solutions at t=0, so a logarithm is forced. The
counting routine is validated on operators with KNOWN log/no-log behaviour first.
"""
import json, os
import sympy as sp

t = sp.symbols('t')

def falling(s, i):
    r = sp.Integer(1)
    for a in range(i):
        r *= (s - a)
    return r

def count_logfree_solutions(op, p, K=8):
    """Number of linearly independent formal solutions y = (t-p)^e * (power series), NO logs,
    within the integer-spaced exponent coset of the smallest exponent. op = [c0..c_n] polynomials.
    Uses exact rational arithmetic; nullity via FLINT (fast, exact) with a sympy fallback."""
    x = sp.symbols('x')
    n = len(op) - 1
    cpoly = [sp.Poly(sp.expand(sp.expand(op[i]).subs(t, p + x)), x) for i in range(n + 1)]
    cij = []
    for i in range(n + 1):
        d = {}
        for (j,), co in cpoly[i].terms():
            if co != 0:                      # skip spurious zero entries from Poly(0,x)
                d[j] = sp.Rational(co)
        cij.append(d)
    # M(x^s)/x^s = sum_d q_d(s) x^d ,  q_d(s) = sum_{i,j: j-i=d} c_{i,j} falling(s,i)
    from collections import defaultdict
    s = sp.symbols('s')
    qd = defaultdict(lambda: sp.Integer(0))
    for i in range(n + 1):
        for j, co in cij[i].items():
            qd[j - i] += co * falling(s, i)
    dmin = min(d for d in qd if sp.expand(qd[d]) != 0)   # lowest power => indicial
    Ipoly = sp.expand(qd[dmin])
    roots = sp.roots(sp.Poly(Ipoly, s))
    exps = []
    for rt, m in roots.items():
        exps += [sp.nsimplify(rt)] * m
    emin = min(exps, key=lambda z: sp.re(sp.N(z)))
    # Frobenius recurrence, LOWER-TRIANGULAR: equation r (coeff of x^{emin+dmin+r}) is
    #   sum_{k<=r} q_{dmin+r-k}(emin+k) a_k = 0 ,  diagonal = q_{dmin}(emin+r) = indicial(emin+r).
    rows = []
    for r in range(K):
        row = []
        for k in range(K):
            d = dmin + r - k
            val = qd[d].subs(s, emin + k) if (k <= r and d in qd) else sp.Integer(0)
            row.append(sp.Rational(sp.expand(val)))
        rows.append(row)
    # nullity = K - rank  (FLINT if available for speed/exactness, else exact sympy)
    try:
        import flint
        Q = flint.fmpq_mat([[flint.fmpq(int(sp.numer(v)), int(sp.denom(v))) for v in r] for r in rows])
        rank = Q.rank()
    except Exception:
        rank = sp.Matrix(rows).rank()
    nullity = K - rank
    return nullity, [str(e) for e in sorted(exps, key=lambda z: sp.re(sp.N(z)))], str(emin)

# ---------- validate on known operators ----------
print("validating count_logfree_solutions on known operators...")
# t^3 D^3 : exponents {0,1,2} at t=0, solutions 1,t,t^2  => 3 log-free (NO log)
op_nolog = [sp.Integer(0), sp.Integer(0), sp.Integer(0), t**3]
d_nolog, e_nolog, _ = count_logfree_solutions(op_nolog, sp.Integer(0))
print("  t^3 D^3  exps=%s  log-free dim=%d  (expect 3, no log)" % (e_nolog, d_nolog))
assert d_nolog == 3, "validation failed (no-log case)"
# theta^3 = t^3 D^3 + 3 t^2 D^2 + t D : exponents {0,0,0}, solutions 1, log, log^2 => 1 log-free
op_log = [sp.Integer(0), t, 3*t**2, t**3]
d_log, e_log, _ = count_logfree_solutions(op_log, sp.Integer(0))
print("  theta^3  exps=%s  log-free dim=%d  (expect 1, logs present)" % (e_log, d_log))
assert d_log == 1, "validation failed (log case)"
# a middle case: theta^2*(theta-1)= exponents {0,0,1}: solutions 1, log, and t => 2 log-free
op_mid = None
print("  [validate] log-detection PASS\n")

# ---------- run on M ----------
HERE = os.path.dirname(os.path.abspath(__file__))
Mc = json.load(open(os.path.join(HERE, "M_coeffs.json")))
M = [sum(int(c) * t**j for j, c in enumerate(Mc["coeffs"][i])) for i in range(4)]

dimM, expsM, eminM = count_logfree_solutions(M, sp.Integer(0))
print("M at t=0 : exponents=%s  emin=%s" % (expsM, eminM))
print("M at t=0 : number of log-free (Laurent) local solutions = %d  (order is 3)" % dimM)
haslog = dimM < 3
print("=> logarithmic solution at t=0 present:", haslog)

print("\n================ NON-LIOUVILLIAN CERTIFICATE ================")
print("1. M is irreducible over Qbar(t)                 : YES (certify_factor.py)")
print("2. M is NOT a symmetric square (excludes SO_3)   : YES (certify_factor.py)")
print("3. M has a genuine logarithmic solution (t=0)    :", haslog)
print("   -> excludes FINITE Galois group (all-algebraic has no logs)")
print("   -> excludes IMPRIMITIVE Galois group (monomial local monodromy is semisimple,")
print("      so an imprimitive operator has no logarithmic solutions)")
print("-------------------------------------------------------------")
if haslog:
    print("=> By Singer-Ulmer, G is neither finite nor imprimitive; G acts irreducibly and")
    print("   non-solvably, so M is NON-LIOUVILLIAN => NO algebraic or elementary closed form.")
    print("   Combined with the SEPARATE not-Sym^2 result (=> no elliptic/K closed form), the")
    print("   hyperkagome LGF has no closed form in algebraic, elementary, OR elliptic functions.")
    print("   (Conditional on M being the minimal annihilator of the LGF; see the module docstring.)")
