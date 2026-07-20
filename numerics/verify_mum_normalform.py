"""
CERTIFY: the local monodromy of M at t=0 is MAXIMALLY UNIPOTENT (MUM), in canonical normal form.

Context.  certify_orthogonal.py part (C) proves the local monodromy at t=0 is a single 3x3
unipotent Jordan block (exponents {-1,-1,0}, all integers; exactly one log-free solution;
maximal log power n=2).  J.-M. Maillard (private communication, July 2026, third message)
observed that the Frobenius basis can moreover be brought to the CANONICAL MUM normal form

    yy0 = y0,
    yy1 = y0*log(t) + f1,
    yy2 = y0*log(t)^2/2 + f1*log(t) + f2,

by the integer recombination yy2 = y2 - (13/30)*y1 of his rescaled basis -- i.e. the SAME
series f1 multiplies log(t) in yy2 as appears in yy1 (the strong Jordan-chain alignment:
under analytic continuation t -> t*e^{2*pi*i}:  yy1 -> yy1 + 2*pi*i*yy0,
yy2 -> yy2 + 2*pi*i*yy1 + (2*pi*i)^2/2 * yy0).  So t=0 is a MUM point in the monodromy sense,
even though the exponents are {-1,-1,0} rather than {0,0,0} (the log partners carry 1/t heads).

This script verifies, in exact rational arithmetic:
  (A) internal consistency of the communicated data: F1 - (13/30)*y0 == f1 on every
      communicated coefficient (including the 1/t heads), so the recombination
      yy2 = y2 - (13/30)*y1 produces the aligned normal form;
  (B) cross-consistency with the first communication: f1 == -(4/15)*(first-message Laurent
      part of y1) and F1 == -(4/15)*(first-message log-coefficient of y2), term by term;
  (C) the DECISIVE ODE identity: M(y0*log t + f1) = 0, i.e. the non-log part
      Cross(y0) + M(f1) vanishes identically through every checkable order, where
      Cross(u) = c1*u/t + c2*(2u'/t - u/t^2) + c3*(3u''/t - 3u'/t^2 + 2u/t^3)
      collects the log-derivative cross terms (D^j log t = (-1)^{j-1}(j-1)! t^{-j}).
      [The same identity with F1 in place of f1 is then automatic, since
       M(F1 - f1) = (13/30) M(y0) = 0.]
Run:  python verify_mum_normalform.py
"""
import json
import os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
Mj = json.load(open(os.path.join(HERE, "M_coeffs.json")))
c = [[F(int(x)) for x in lst] for lst in Mj["coeffs"]]          # c0..c3 ascending

nu = [F(a, b) for a, b in json.load(open(os.path.join(HERE, "nu.json")))]
y0 = [F(n + 1) * nu[n + 1] / 2 for n in range(len(nu) - 1)]     # analytic solution, ~109 terms

# --- data communicated by Maillard (message 3), coefficients of t^{-1}..t^{10} ---
f1 = [F(4, 15), F(8, 3), F(251, 10), F(1092, 5), F(18457, 10), F(231599, 15),
      F(3872419, 30), F(16250528, 15), F(1920237593, 210), F(2717207602, 35),
      F(417545285593, 630), F(358201935913, 63)]
F1 = [F(4, 15), F(31, 10), F(883, 30), F(2561, 10), F(64783, 30), F(540613, 30),
      F(4510537, 30), F(37799713, 30), F(2230870097, 210), F(6308648291, 70),
      F(484447525543, 630), F(4154283780859, 630)]
# --- message 1 data (see _agents/maillard_email_technical.md), same index range ---
L1_old = [F(-1), F(-10), F(-753, 8), F(-819), F(-55371, 8), F(-231599, 4), F(-3872419, 8),
          F(-4062632), F(-1920237593, 56), F(-4075811403, 14), F(-417545285593, 168),
          F(-1791009679565, 84)]
L2_old = [F(-1), F(-93, 8), F(-883, 8), F(-7683, 8), F(-64783, 8), F(-540613, 8),
          F(-4510537, 8), F(-37799713, 8), F(-2230870097, 56), F(-18925944873, 56),
          F(-484447525543, 168), F(-4154283780859, 168)]

# (A) recombination alignment: F1 - (13/30)*y0 == f1  (1/t heads carry no y0 term)
okA = (F1[0] == f1[0]) and all(F1[k + 1] - F(13, 30) * y0[k] == f1[k + 1] for k in range(11))
print("(A) F1 - (13/30)*y0 == f1 on all communicated terms:", "PASS" if okA else "FAIL")

# (B) rescaling vs message 1: f1 == -(4/15)*L1_old,  F1 == -(4/15)*L2_old
okB = all(f1[k] == F(-4, 15) * L1_old[k] for k in range(12)) and \
      all(F1[k] == F(-4, 15) * L2_old[k] for k in range(12))
print("(B) message-3 series == -(4/15) * message-1 series, all terms:", "PASS" if okB else "FAIL")

# (C) M(y0*log t + f1) = 0: the non-log part  Cross(y0) + M(f1)  vanishes.
def deriv(a, off):
    """(coeffs, offset) Laurent derivative."""
    return [a[i] * F(off + i) for i in range(len(a))], off - 1

def shift(a, off, k):
    return a, off + k

def polymul_ser(pol, ser, off):
    r = [F(0)] * (len(pol) + len(ser) - 1)
    for i, pi in enumerate(pol):
        if pi:
            for j, sj in enumerate(ser):
                if sj:
                    r[i + j] += pi * sj
    return r, off

acc = {}
def deposit(ser, off, limit):
    for i, v in enumerate(ser):
        k = i + off
        if v and k <= limit:
            acc[k] = acc.get(k, F(0)) + v

# validity limit: f1 known through t^10  =>  c0*f1 valid to 10; c1*f1' to 9; c2 (val 1)*f1'' to 9;
# c3 (val 2)*f1''' to 9.  Cross(y0) is valid far beyond.  Check through t^9.
LIMIT = 9
u = y0; uo = 0
u1, u1o = deriv(u, uo)
u2, u2o = deriv(u1, u1o)
# Cross(y0) = c1*u/t + c2*(2u'/t - u/t^2) + c3*(3u''/t - 3u'/t^2 + 2u/t^3)
for pol, ser, off in [
    (c[1], u, uo - 1),
    (c[2], [2 * x for x in u1], u1o - 1), (c[2], [-x for x in u], uo - 2),
    (c[3], [3 * x for x in u2], u2o - 1), (c[3], [-3 * x for x in u1], u1o - 2),
    (c[3], [2 * x for x in u], uo - 3),
]:
    s, o = polymul_ser(pol, ser, off)
    deposit(s, o, LIMIT)
# + M(f1):  f1 has offset -1
g = f1; go = -1
g1, g1o = deriv(g, go)
g2, g2o = deriv(g1, g1o)
g3, g3o = deriv(g2, g2o)
for pol, ser, off in [(c[0], g, go), (c[1], g1, g1o), (c[2], g2, g2o), (c[3], g3, g3o)]:
    s, o = polymul_ser(pol, ser, off)
    deposit(s, o, LIMIT)
bad = sorted(k for k, v in acc.items() if v != 0)
okC = not bad
lo = min(acc.keys()) if acc else None
print("(C) M(y0*log t + f1) == 0: non-log part vanishes for all orders t^%s..t^%d:" % (lo, LIMIT),
      "PASS" if okC else ("FAIL at %s" % bad[:6]))

print()
if okA and okB and okC:
    print("================ MUM NORMAL-FORM CERTIFICATE ================")
    print("The point t=0 is a MUM point of M (maximally unipotent local monodromy, single 3x3")
    print("Jordan block -- certify_orthogonal.py (C)), and the Frobenius basis takes the")
    print("canonical MUM normal form  {y0,  y0*log t + f1,  y0*log^2(t)/2 + f1*log t + f2}")
    print("after the integer recombination yy2 = y2 - (13/30)*y1")
    print("[J.-M. Maillard, private communication, July 2026; verified here exactly].")
    print("The exponents are {-1,-1,0} (not {0,0,0}): MUM in the monodromy sense, with 1/t")
    print("Laurent heads on the log partners -- t=0 is E=infinity in the energy variable.")
else:
    raise SystemExit(1)
