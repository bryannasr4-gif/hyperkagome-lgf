"""Independent numerical confirmation of the exact special value Re G(1) = 1/9, by direct
Brillouin-zone averaging of the resolvent over the 12x12 hyperkagome Bloch Hamiltonian.

Re G(1) = flat-band pole (1/3)/(1-(-2)) = 1/9  PLUS a dispersive principal value that
vanishes by the exact E=1 reflection symmetry.  The cancellation is exact only under a
symmetric quadrature: a naive Gamma-inclusive grid samples the E=1 van Hove states unevenly
and leaves a spurious few-percent offset, so we use an off-Gamma (half-shifted) grid.
"""
import pickle, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sites, bonds = pickle.load(open(os.path.join(HERE, "lattice.pkl"), "rb"))
B = [(int(i), int(j), tuple(int(round(x)) for x in d)) for (i, j, d) in bonds]
nsite = 12

def eigs(N, shift):
    ks = (np.arange(N) + shift) / N * 2 * np.pi
    KX, KY, KZ = np.meshgrid(ks, ks, ks, indexing='ij')
    kx, ky, kz = KX.ravel(), KY.ravel(), KZ.ravel()
    H = np.zeros((kx.size, nsite, nsite), dtype=complex)
    for (a, b, d) in B:
        H[:, a, b] += np.exp(1j * (kx*d[0] + ky*d[1] + kz*d[2]))
    return np.linalg.eigvalsh((H + np.conjugate(np.transpose(H, (0, 2, 1)))) / 2)

def ReG(z, N, shift):
    ev = eigs(N, shift)
    return float(np.mean(np.real(1.0 / (z - ev))))

z = 1.0 + 1e-9j
print("Re G(1) = 1/9 = %.10f (exact)\n" % (1/9))
print(" grid   off-Gamma (shift 0.5)     naive Gamma-incl (shift 0.0)")
for N in [24, 40, 60, 84]:
    print("  %2d^3   %.6f                 %.6f" % (N, ReG(z, N, 0.5), ReG(z, N, 0.0)))
val = ReG(z, 84, 0.5)
print("\noff-Gamma Re G(1) -> %.6f   (target 1/9 = %.6f, |err| = %.2e)"
      % (val, 1/9, abs(val - 1/9)))
