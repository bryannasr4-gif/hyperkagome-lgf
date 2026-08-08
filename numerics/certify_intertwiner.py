"""
CERTIFY: why the intertwiner between M and Sym^2(V_2) is invisible over Q(t), and
         exactly which conjugation makes it appear.

MOTIVATION (J.-M. Maillard, private communication, August 2026).  Searching directly for a
homomorphism between the hyperkagome operator M and the symmetric square of V_2 in its
PROJECTIVE NORMAL FORM

        V_2 : w'' + Q_V w = 0,       Ntilde := Sym^2(V_2) = D^3 + 4 Q_V D + 2 Q_V',

returns NOTHING in either direction -- Maple's Homomorphisms(symmetric_power(V2,2), M) and
Homomorphisms(M, symmetric_power(V2,2)) are both empty.  That is not a limitation of the
implementation.  This certificate proves it is FORCED, identifies the obstruction as the
determinant (Wronskian) character of M, and shows that conjugating M by the twist unit
v = sqrt((1-4t)(1-5t)(1-9t)) removes exactly that obstruction and nothing else.

WHAT IS PROVED HERE (all exact, in Q(t); no floating point anywhere)

  [1] The Wronskian log-derivative of M, a := -c2/c3, has HALF-INTEGER residues at
      t = 1/4, 1/5, 1/9 (and integer residues elsewhere).  Hence a is not f'/f for any
      f in Q(t)^*, so the determinant module Lambda^3(M) is NOT trivial over Q(t).

  [2] Ntilde has c2 = 0 identically, so Lambda^3(Ntilde) IS trivial.

  [3] Therefore Lambda^3(M) and Lambda^3(Ntilde) are non-isomorphic, so M and Ntilde are
      non-isomorphic differential modules over Q(t).  Both being irreducible (an input:
      certify_factor.py for M; for Ntilde, Sym^2 of the uniformizing operator of
      X(Gamma_0(30)^+), whose projective monodromy is Zariski-dense in PSL(2,C)), every
      nonzero homomorphism between them would be an isomorphism.  There is none, in either
      direction.  This is the empty Homomorphisms output, explained.

  [4] M_v := v o M o v^{-1} has coefficients in Q(t) -- verified coefficient by coefficient,
      not assumed -- and its Wronskian log-derivative IS a log-derivative of an explicit
      rational function, exhibited here.  Why the character dies needs the right premise, and
      it is NOT "the character has order two" alone: characters multiply, so conjugation sends
      chi to chi * chi_v^3 = chi * chi_v, which is trivial precisely when chi = chi_v.  That
      equality is the certified content of certify_orthogonal.py stage (D): the determinant
      character of M IS the quadratic character of v^2 = (1-4t)(1-5t)(1-9t).  Given it,
      chi * chi_v = chi_v^2 = 1.  The explicit witness below proves the trivialization outright
      and does not depend on that bookkeeping.

  [5] Over M_v the intertwiner is rational and has ORDER ONE: with T = rho_0 + rho_1 D of
      Eq. (15) of the paper, the right-division remainder of M_v o T by Ntilde is exactly 0,
      so T maps Sol(Ntilde) into Sol(M_v).  (Independently recomputed here; certify_y0.py
      reaches the same identity by a different route.)

  [6] The obstruction is an order-two character, so it dies in every even tensor
      construction: 2a and 4a have integer residues where a does not.  That is why the
      determinant obstruction is absent at the symmetric-square level -- Sym^2(M) against
      Sym^4(V_2) -- where the direct search at the level of M is blocked.  (The homomorphy
      at that level was reported by Maillard's Maple session and is NOT certified here;
      only the disappearance of the obstruction is.)

NEGATIVE CONTROLS (a check that has never been shown to fail is not a check)
  - the log-derivative predicate REJECTS a function with a half-integer residue, REJECTS a
    function with a double pole, and REJECTS one carrying a polynomial part;
  - right division REJECTS a non-factor;
  - perturbing rho_1 by t breaks the intertwiner identity of [5];
  - positive witness controls (a rational function with an irreducible degree-7 factor
    included) are correctly ACCEPTED, so [1] is not an artifact of the routine always
    saying "no".

Runs on system Python (sympy only).  Writes CERTIFICATE_intertwiner.txt.
"""
import json
import os
import sys

import sympy as sp

t = sp.Symbol('t')
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = []
FAILS = []


def out(s=""):
    print(s)
    OUT.append(s)


def check(name, ok, detail=""):
    out("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("   " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)
    return ok


# ----------------------------------------------------------------- operator algebra
# An operator is a list [a0, a1, ...] meaning sum_i a_i D^i, coefficients in Q(t).

def trim(a):
    a = list(a)
    while len(a) > 1 and sp.simplify(a[-1]) == 0:
        a.pop()
    return a


def opadd(a, b):
    n = max(len(a), len(b))
    return trim([sp.cancel((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0))
                 for i in range(n)])


def opscale(a, c):
    return trim([sp.cancel(c * x) for x in a])


def opmul(a, b):
    """Composition a o b in Q(t)<D>, D f = f D + f'."""
    res = [sp.Integer(0)]
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        # D^i applied on the left of b: Leibniz
        cur = list(b)
        for _ in range(i):
            nxt = [sp.Integer(0)] * (len(cur) + 1)
            for k, ck in enumerate(cur):
                nxt[k + 1] += ck                      # D * (c D^k) = c D^{k+1} + c' D^k
                nxt[k] += sp.diff(ck, t)
            cur = [sp.cancel(x) for x in nxt]
        res = opadd(res, opscale(cur, ai))
    return trim(res)


def rdiv(a, b):
    """Right division a = Q o b + R with order(R) < order(b).  Returns (Q, R)."""
    a = [sp.cancel(x) for x in list(a)]
    db = len(b) - 1
    lb = b[-1]
    Q = [sp.Integer(0)]
    while len(trim(a)) - 1 >= db and trim(a) != [sp.Integer(0)]:
        a = trim(a)
        da = len(a) - 1
        if da < db:
            break
        c = sp.cancel(a[-1] / lb)
        shift = [sp.Integer(0)] * (da - db) + [c]
        Q = opadd(Q, shift)
        a = opadd(a, opscale(opmul(shift, b), -1))
    return trim(Q), trim(a)


# ---------------------------------------------------------- residues and witnesses
#
# The two directions are deliberately asymmetric, because they need different things.
#
#   NEGATIVE ("g is not a log-derivative"):  a SINGLE exact residue suffices and is a
#   complete proof.  Every residue of f'/f is an integer, for any f in Qbar(t)^* (near a
#   zero or pole of order n, f'/f = n/(t-a) + O(1)).  So one non-integer residue at one
#   rational point rules out every f, over Qbar(t) and a fortiori over Q(t).  No
#   enumeration of the other poles is involved, so nothing can be silently missed.
#
#   POSITIVE ("g is a log-derivative"):  an explicit witness f is constructed and the
#   identity f'/f == g is then verified exactly.  If the construction fails we report
#   "no witness", which is NOT a claim of impossibility.
#
# An earlier version used sp.roots() to enumerate poles.  That silently skipped the poles
# on the degree-7 apparent locus p_7, which sp.roots cannot split, and the reconstruction
# step is what exposed it.  The rewrite below never enumerates roots.

def residue_at(g, pt):
    """Exact residue of g in Q(t) at the rational point pt."""
    return sp.nsimplify(sp.simplify(sp.residue(sp.cancel(g), t, sp.nsimplify(pt))))


def log_derivative_witness(g):
    """Try to write g = f'/f with f in Q(t)^*, by partial fractions over Q (never over
    Qbar, so the irreducible factors of the denominator stay intact).  Returns (f, True)
    only after verifying f'/f - g == 0 exactly; otherwise (reason, False)."""
    g = sp.cancel(sp.together(g))
    num, den = sp.fraction(g)
    dpoly = sp.Poly(den, t)
    if sp.degree(sp.Poly(num, t)) >= dpoly.degree():
        quo, _ = sp.div(sp.Poly(num, t), dpoly)
        if quo.as_expr() != 0:
            return "nonzero polynomial part %s" % quo.as_expr(), False
    facs = sp.factor_list(dpoly.as_expr())[1]
    f = sp.Integer(1)
    Np = sp.Poly(num, t, domain='QQ')
    for P, mult in facs:
        if mult != 1:
            return "denominator factor %s has multiplicity %d" % (sp.factor(P), mult), False
        # The multiplicity n attached to the irreducible P is read off in the residue field
        # Q[t]/(P):  g = n P'/P + (regular), so  N == n P' (D/P)  mod P, i.e.
        #            n == N * (D/P)^{-1} * (P')^{-1}   in Q[t]/(P).
        # Nothing here enumerates the roots of P, so a high-degree irreducible factor such
        # as the apparent locus p_7 is handled exactly like a linear one.
        Pp = sp.Poly(P, t, domain='QQ')
        Dip = sp.Poly(sp.cancel(dpoly.as_expr() / P), t, domain='QQ')
        dPp = sp.Poly(sp.diff(P, t), t, domain='QQ')
        try:
            n_poly = (Np * sp.invert(Dip, Pp) * sp.invert(dPp, Pp)) % Pp
        except sp.polys.polyerrors.NotInvertible:
            return "non-invertible cofactor on %s" % sp.factor(P), False
        if n_poly.degree() > 0:
            return "multiplicity on %s is not constant" % sp.factor(P), False
        n_val = sp.Rational(n_poly.as_expr())
        if not n_val.is_Integer:
            return "non-integer multiplicity %s on factor %s" % (n_val, sp.factor(P)), False
        f *= P ** int(n_val)
    if sp.simplify(sp.cancel(sp.diff(f, t) / f - g)) != 0:
        return "witness did not verify", False
    return sp.factor(f), True


out("=" * 78)
out("certify_intertwiner.py -- the determinant-character obstruction to a rational")
out("                          intertwiner between M and Sym^2(V_2), and its removal")
out("=" * 78)

# =========================================================== Stage 0: primitives
out("\n--- Stage 0: validation of the primitives on objects of KNOWN structure ------")

# operator composition against a hand-computed case
lhs = opmul([sp.Integer(0), sp.Integer(1)], [t, sp.Integer(1)])       # D o (t + D)
check("D o (t + D) = 1 + t D + D^2",
      [sp.simplify(x) for x in lhs] == [sp.Integer(1), t, sp.Integer(1)], str(lhs))

# right division recovers a factor, and rejects a non-factor
B = [sp.Integer(-2), sp.Integer(1)]                                    # D - 2
A = opmul([t, sp.Integer(1)], B)                                       # (t + D) o (D - 2)
_, R = rdiv(A, B)
check("right division of a genuine multiple leaves remainder 0", trim(R) == [sp.Integer(0)])
_, R2 = rdiv(opadd(A, [sp.Integer(1)]), B)
check("NEGATIVE CONTROL: right division of a NON-multiple leaves remainder != 0",
      trim(R2) != [sp.Integer(0)], "remainder %s" % R2)

# the witness constructor, on a control whose answer is known
f_ctl = (t - 1) ** 2 / t ** 3
w_ctl, ok = log_derivative_witness(sp.cancel(sp.diff(f_ctl, t) / f_ctl))
check("witness constructor RECOVERS f = (t-1)^2/t^3 from f'/f", ok, str(w_ctl))

# and on a control carrying an IRREDUCIBLE high-degree factor -- the case that broke an
# earlier root-enumerating version of this routine (the degree-7 apparent locus p_7).
p7 = (101025 * t ** 7 - 369600 * t ** 6 + 455798 * t ** 5 - 290956 * t ** 4
      + 93657 * t ** 3 - 17580 * t ** 2 + 1800 * t - 64)
check("p_7 is irreducible over Q (so its roots cannot be enumerated rationally)",
      len(sp.factor_list(p7)[1]) == 1 and sp.factor_list(p7)[1][0][1] == 1)
f_p7 = p7 ** 3 / t ** 2
w_p7, ok = log_derivative_witness(sp.cancel(sp.diff(f_p7, t) / f_p7))
check("witness constructor RECOVERS f = p_7^3/t^2, poles on p_7 included", ok, str(w_p7)[:70])

w_bad, ok = log_derivative_witness(sp.Rational(1, 2) / (t - 1))
check("NEGATIVE CONTROL: no witness for a half-integer multiplicity", not ok, str(w_bad))
w_bad, ok = log_derivative_witness(1 / (t - 1) ** 2)
check("NEGATIVE CONTROL: no witness for a double pole", not ok, str(w_bad))
w_bad, ok = log_derivative_witness(1 / (t - 1) + t)
check("NEGATIVE CONTROL: no witness for a nonzero polynomial part", not ok, str(w_bad))

# the residue routine, on which every NEGATIVE claim below rests
check("residue routine: res of 3/(t-1) at t=1 is 3", residue_at(3 / (t - 1), 1) == 3)
check("residue routine: res of -5/(2(4t-1)) at t=1/4 is -5/8",
      residue_at(-sp.Rational(5, 2) / (4 * t - 1), sp.Rational(1, 4)) == sp.Rational(-5, 8))
check("residue routine: res at a regular point is 0", residue_at(1 / (t - 1), 5) == 0)

# =========================================================== Stage 1: M's character
out("\n--- Stage 1: the determinant character of M ----------------------------------")
Md = json.load(open(os.path.join(HERE, "M_coeffs.json")))
assert Md["order"] == 3, "M_coeffs.json is not order 3"
Mcoef = [sum(sp.Integer(c) * t ** k for k, c in enumerate(col)) for col in Md["coeffs"]]
check("M read from M_coeffs.json has order 3 and nonzero leading coefficient",
      len(Mcoef) == 4 and sp.simplify(Mcoef[3]) != 0)

a = sp.cancel(-Mcoef[2] / Mcoef[3])          # W'/W for the Wronskian W of M
SING = [sp.Integer(0), sp.Rational(1, 9), sp.Rational(1, 5), sp.Rational(1, 4), sp.Integer(1)]
res_map = {r: residue_at(a, r) for r in SING}
out("    Wronskian log-derivative a = -c2/c3 ; residues at the finite singular points:")
for r in SING:
    out("       t = %-10s residue %s" % (r, res_map[r]))

half = [sp.Rational(1, 4), sp.Rational(1, 5), sp.Rational(1, 9)]
check("the residues at t = 1/4, 1/5, 1/9 are all HALF-INTEGERS (not integers)",
      all(sp.Rational(res_map[r]).q == 2 for r in half),
      ", ".join("res(%s)=%s" % (r, res_map[r]) for r in half))
check("the residues at t = 0 and t = 1 ARE integers (so the obstruction is localized)",
      all(sp.Rational(res_map[r]).is_Integer for r in (sp.Integer(0), sp.Integer(1))))
# One non-integer residue is a COMPLETE proof: every residue of f'/f is an integer.
ok_a = sp.Rational(res_map[sp.Rational(1, 4)]).is_Integer
check("hence a is NOT f'/f for ANY f in Qbar(t)* => Lambda^3(M) is non-trivial",
      not ok_a, "res at t=1/4 is %s, not an integer" % res_map[sp.Rational(1, 4)])
w_a, ok_wit = log_derivative_witness(a)
check("cross-check: the witness constructor also finds no witness for a", not ok_wit, str(w_a))

# =========================================================== Stage 2: Ntilde
out("\n--- Stage 2: Sym^2(V_2) in projective normal form is unimodular --------------")
Vd = json.load(open(os.path.join(HERE, "V2_data.json")))
QV = sp.cancel(sp.sympify(Vd["Q"]["num"]) / sp.sympify(Vd["Q"]["den"]))
Ntil = [sp.cancel(2 * sp.diff(QV, t)), sp.cancel(4 * QV), sp.Integer(0), sp.Integer(1)]
check("Ntilde = D^3 + 4 Q_V D + 2 Q_V' has a vanishing D^2 coefficient",
      sp.simplify(Ntil[2]) == 0)
check("so its Wronskian log-derivative is 0, a log-derivative (f = 1): Lambda^3 trivial",
      log_derivative_witness(sp.Integer(0))[1])
# sanity: Q_V really is the projective normal form of the stored V2 = D^2 + p D + q
p_ = sp.cancel(sp.sympify(Vd["p"]["num"]) / sp.sympify(Vd["p"]["den"]))
q_ = sp.cancel(sp.sympify(Vd["q"]["num"]) / sp.sympify(Vd["q"]["den"]))
check("cross-check: Q_V == q - p'/2 - p^2/4 from the stored V_2",
      sp.simplify(QV - (q_ - sp.diff(p_, t) / 2 - p_ ** 2 / 4)) == 0)

# =========================================================== Stage 3: the conclusion
out("\n--- Stage 3: no rational intertwiner, in either direction --------------------")
check("Lambda^3(M) non-trivial and Lambda^3(Ntilde) trivial => M is NOT isomorphic to Ntilde",
      (not ok_a) and sp.simplify(Ntil[2]) == 0)
out("    Both operators are irreducible (M: certify_factor.py; Ntilde: Sym^2 of the")
out("    uniformizing operator of X(Gamma_0(30)^+), projective monodromy Zariski-dense in")
out("    PSL(2,C)), so a nonzero homomorphism would be an isomorphism.  Hence")
out("    Hom(M, Ntilde) = Hom(Ntilde, M) = 0 over Q(t): the empty Maple output is FORCED.")

# =========================================================== Stage 4: the conjugation
out("\n--- Stage 4: conjugating by v = sqrt((1-4t)(1-5t)(1-9t)) kills the character --")
w = sp.cancel(sp.diff(sp.log((1 - 4 * t) * (1 - 5 * t) * (1 - 9 * t)), t) / 2)   # v'/v
check("v'/v is RATIONAL (so the conjugation stays inside Q(t)<D>)",
      sp.simplify(sp.together(w).is_rational_function(t)) is not False, str(sp.factor(w)))

# M_v = v o M o v^{-1} = sum_i c_i (D - w)^i
Dm = [sp.cancel(-w), sp.Integer(1)]
powers = [[sp.Integer(1)]]
for _ in range(3):
    powers.append(opmul(powers[-1], Dm))
Mv = [sp.Integer(0)]
for i in range(4):
    Mv = opadd(Mv, opscale(powers[i], Mcoef[i]))
Mv = trim(Mv)
check("M_v := v o M o v^{-1} has order 3", len(Mv) - 1 == 3)
check("every coefficient of M_v lies in Q(t) (verified one by one, not assumed)",
      all(sp.together(c).is_rational_function(t) for c in Mv))

av = sp.cancel(-Mv[2] / Mv[3])
for r in half:
    check("the residue of M_v's character at t = %s is now an INTEGER" % r,
          sp.Rational(residue_at(av, r)).is_Integer, "residue %s" % residue_at(av, r))
f_v, ok_v = log_derivative_witness(av)
check("an explicit rational witness f with f'/f = -c2/c3 exists => character KILLED",
      ok_v, "f = %s" % str(f_v)[:90])
check("and it differs from a by exactly 3 v'/v, the conjugation shift",
      sp.simplify(av - (a + 3 * w)) == 0)

# =========================================================== Stage 5: the intertwiner
out("\n--- Stage 5: over M_v the intertwiner is rational and of ORDER ONE -----------")
rho1 = sp.cancel((15 * t ** 2 + 17 * t - 8) / (30 * t * (t - 1)))
rho0 = sp.cancel(-(4050 * t ** 6 + 2445 * t ** 5 - 11436 * t ** 4 + 8000 * t ** 3
                   - 2130 * t ** 2 + 231 * t - 8)
                 / (30 * t ** 2 * (t - 1) ** 2 * (4 * t - 1) * (5 * t - 1) * (9 * t - 1)))
T = [rho0, rho1]
check("T = rho_0 + rho_1 D has order ONE", len(trim(T)) - 1 == 1)
_, Rem = rdiv(opmul(Mv, T), Ntil)
check("right-division remainder of M_v o T by Ntilde is exactly 0 "
      "=> T : Sol(Ntilde) -> Sol(M_v)", trim(Rem) == [sp.Integer(0)], "remainder %s" % Rem)
_, RemBad = rdiv(opmul(Mv, [rho0, sp.cancel(rho1 + t)]), Ntil)
check("NEGATIVE CONTROL: perturbing rho_1 by t breaks that identity",
      trim(RemBad) != [sp.Integer(0)])
_, RemNoConj = rdiv(opmul(Mcoef, T), Ntil)
check("NEGATIVE CONTROL: the SAME T against the UNCONJUGATED M leaves a nonzero remainder",
      trim(RemNoConj) != [sp.Integer(0)])

# =========================================================== Stage 6: even powers
out("\n--- Stage 6: why the symmetric-square level sees what M does not -------------")
out("    For a rank-three module, det Sym^2 = (det)^4 and det Sym^2(M) is therefore the")
out("    fourth power of the character; the character has order two, so every even tensor")
out("    power of it is trivial.  Concretely:")
for k in (1, 2, 4):
    r14 = residue_at(sp.cancel(k * a), sp.Rational(1, 4))
    integral = sp.Rational(r14).is_Integer
    check("%d * a has %s residue at t = 1/4  (%s)"
          % (k, "an INTEGER" if integral else "a non-integer", r14),
          integral == (k % 2 == 0))
wit2, ok2 = log_derivative_witness(sp.cancel(2 * a))
check("and 2a really is a log-derivative: explicit witness verified",
      ok2, "f = %s" % str(wit2)[:70])

# ============================== Stage 7: consistency with the rational bridge operator
out("\n--- Stage 7: consistency with the rational P of the bridge identity ----------")
out("    The paper also exhibits a RATIONAL order-two P over Q(t) with P(y_0) = f_0^2,")
out("    M(y_0) = 0 and V_2(f_0) = 0.  Taken carelessly that contradicts Stage 3.  It does")
out("    not, because there V_2 is the stored NON-unimodular form D^2 + p D + q, whose")
out("    symmetric square has Wronskian log-derivative -3p, not 0.  Checked here:")
a_stored = sp.cancel(-3 * p_)          # -c2/c3 of Sym^2(D^2 + p D + q)
for r in half:
    rr = residue_at(a_stored, r)
    check("res of the STORED Sym^2(V_2) character at t = %s is %s, also a half-integer"
          % (r, rr), sp.Rational(rr).q == 2)
check("so it does NOT agree with the unimodular case (which has residue 0 everywhere)",
      any(residue_at(a_stored, r) != 0 for r in half))
f_cons, ok_cons = log_derivative_witness(sp.cancel(a - a_stored))
check("M and the STORED Sym^2(V_2) have ISOMORPHIC determinants: explicit rational witness",
      ok_cons, "f = %s" % str(f_cons)[:80])

# Equal determinants only remove the obstruction.  Exhibit the intertwiner itself: with
# N the stored order-three operator and P = v0 + v1 D + v2 D^2 the conic/bridge operator,
# P is a module homomorphism iff N o P vanishes on Sol(M), i.e. iff the right-division
# remainder of N o P by M is exactly zero.  (certify_bridge.py verifies only the
# single-solution identity f_0^2 = P(y_0); this is the stronger operator statement, and it
# is what justifies calling P an intertwiner.)
Nst = [sp.cancel(sp.sympify(Vd[k]["num"]) / sp.sympify(Vd[k]["den"]))
       for k in ("B0", "B1", "B2")] + [sp.Integer(1)]
Pop = [sp.cancel(sp.sympify(Vd[k]["num"]) / sp.sympify(Vd[k]["den"]))
       for k in ("v0", "v1", "v2")]
check("the stored order-three operator N really is Sym^2(V_2): B0 == 4pq + 2q'",
      sp.simplify(Nst[0] - (4 * p_ * q_ + 2 * sp.diff(q_, t))) == 0)
check("... and B1 == 2p^2 + p' + 4q",
      sp.simplify(Nst[1] - (2 * p_ ** 2 + sp.diff(p_, t) + 4 * q_)) == 0)
check("... and B2 == 3p", sp.simplify(Nst[2] - 3 * p_) == 0)
_, RemP = rdiv(opmul(Nst, Pop), Mcoef)
check("P is a genuine INTERTWINER: remainder of (N o P) mod M is exactly 0",
      trim(RemP) == [sp.Integer(0)], "order(P) = %d" % (len(trim(Pop)) - 1))
_, RemPbad = rdiv(opmul(Nst, [Pop[0], Pop[1], sp.cancel(Pop[2] + t)]), Mcoef)
check("NEGATIVE CONTROL: perturbing the top coefficient of P breaks it",
      trim(RemPbad) != [sp.Integer(0)])
_, RemPproj = rdiv(opmul(Ntil, Pop), Mcoef)
check("NEGATIVE CONTROL: the same P against the PROJECTIVE Sym^2(V_2) does NOT intertwine",
      trim(RemPproj) != [sp.Integer(0)])
out("    Hence a rational intertwiner is possible in the stored normalization and impossible")
out("    in the projective one.  Both of the paper's statements stand; the distinction is")
out("    exactly the normalization of V_2, which is why the passage needed spelling out.")

# =========================================================== result
out("\n" + "=" * 78)
if FAILS:
    out("RESULT: %d CHECK(S) FAILED" % len(FAILS))
    for f in FAILS:
        out("   - %s" % f)
else:
    out("RESULT: ALL CHECKS PASS.")
    out("The determinant character of M -- the quadratic character of")
    out("v^2 = (1-4t)(1-5t)(1-9t) -- is the ONLY obstruction to a rational intertwiner")
    out("with Sym^2(V_2) in projective normal form.  It makes Hom(M, Sym^2(V_2)) = 0 over")
    out("Q(t) in both directions, it is removed by the single conjugation M_v = v M v^{-1},")
    out("and over M_v the intertwiner is the ORDER-ONE operator T = rho_0 + rho_1 d/dt.")
    out("Being of order two, the character is invisible at the symmetric-square level,")
    out("which is why the determinant obstruction is absent between Sym^2(M) and")
    out("Sym^4(V_2) (their homomorphy was reported by J.-M. Maillard, private")
    out("communication, and is NOT certified by this script) while M and Sym^2(V_2)")
    out("are provably not homomorphic.")
out("=" * 78)

open(os.path.join(HERE, "CERTIFICATE_intertwiner.txt"), "w").write("\n".join(OUT) + "\n")
sys.exit(1 if FAILS else 0)
