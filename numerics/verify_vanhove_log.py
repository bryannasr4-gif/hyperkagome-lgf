"""Numerical support for the van Hove statement at E = 2 (equivalently E = 0): the
Lorentzian-broadened density of states grows linearly in log(1/eta), as a logarithmic
divergence requires, while control points saturate.

Background: t = (E-1)^{-2} maps BOTH E = 0 and E = 2 to the singular point t = 1 of the
operator M, whose Frobenius exponents there are {-1, -1/2, +1/2}. The integer-separated
pair (-1/2, +1/2) permits a logarithmic local solution; this script checks that the
spectrum actually behaves that way, and that E = 0 and E = 2 agree, as the exact E = 1
reflection symmetry requires.

Method: rho_eta(E) = -(1/pi) Im mean_k tr (E + i eta - H(k))^{-1} is computed on an
off-Gamma 64^3 Brillouin-zone grid for a ladder of halving eta. For a log divergence the
increment rho_{eta/2} - rho_eta tends to a positive constant (A log 2); for a finite DOS
it tends to zero. Criteria (with margin; smallest eta excluded as it saturates on the
finite grid):
  (1) at E = 2 every increment over eta = 0.1 -> 0.0125 exceeds 0.03,
  (2) at the smooth control E = 3 every increment is below 0.025,
  (3) rho(E=0) and rho(E=2) agree to 5e-3 at every eta (same point t = 1).

This is a numerical consistency check, not a proof; the paper states only that the local
exponents permit the logarithm and that this broadening test is consistent with it.
"""
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sites, bonds = pickle.load(open(os.path.join(HERE, "lattice.pkl"), "rb"))
B = [(int(i), int(j), tuple(int(round(x)) for x in d)) for (i, j, d) in bonds]
nsite = 12

N = 64
CHUNK = 4096
ETAS = [0.2, 0.1, 0.05, 0.025, 0.0125]
PROBES = [0.0, 2.0, 3.0]

ks = (np.arange(N) + 0.5) / N * 2 * np.pi
KX, KY, KZ = np.meshgrid(ks, ks, ks, indexing='ij')
kx, ky, kz = KX.ravel(), KY.ravel(), KZ.ravel()
npts = kx.size

acc = {(E, eta): 0.0 for E in PROBES for eta in ETAS}
for s in range(0, npts, CHUNK):
    e = slice(s, min(s + CHUNK, npts))
    H = np.zeros((kx[e].size, nsite, nsite), dtype=complex)
    for (a, b, d) in B:
        H[:, a, b] += np.exp(1j * (kx[e] * d[0] + ky[e] * d[1] + kz[e] * d[2]))
    ev = np.linalg.eigvalsh((H + np.conjugate(np.transpose(H, (0, 2, 1)))) / 2).ravel()
    for E in PROBES:
        d2 = (E - ev) ** 2
        for eta in ETAS:
            acc[(E, eta)] += float(np.sum(eta / (d2 + eta * eta)))

ntot = npts * nsite
rho = {k: v / ntot / np.pi for k, v in acc.items()}

print("Lorentzian-broadened DOS (64^3 off-Gamma grid, per site):")
print("  eta      " + "".join("     E=%-6.1f" % E for E in PROBES))
for eta in ETAS:
    print("  %-8.4f" % eta + "".join("  %10.5f" % rho[(E, eta)] for E in PROBES))

ok = True

def check(label, cond):
    global ok
    print("  %-64s %s" % (label, "PASS" if cond else "FAIL"))
    ok = ok and cond

print()
inc2 = [rho[(2.0, ETAS[i + 1])] - rho[(2.0, ETAS[i])] for i in range(len(ETAS) - 1)]
inc3 = [rho[(3.0, ETAS[i + 1])] - rho[(3.0, ETAS[i])] for i in range(len(ETAS) - 1)]
print("  increments per halving  E=2:", ["%.4f" % v for v in inc2])
print("  increments per halving  E=3:", ["%.4f" % v for v in inc3])
check("(1) E=2 increments all exceed 0.03 (log-divergent growth)",
      all(v > 0.03 for v in inc2))
check("(2) E=3 (smooth control) increments all below 0.025",
      all(v < 0.025 for v in inc3))
check("(3) rho(E=0) == rho(E=2) to 5e-3 at every eta (same point t=1)",
      all(abs(rho[(0.0, eta)] - rho[(2.0, eta)]) < 5e-3 for eta in ETAS))

print()
print("ALL CHECKS PASS" if ok else "FAILURE")
sys.exit(0 if ok else 1)
