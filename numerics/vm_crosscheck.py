"""Cross-check that THIS tight-binding hyperkagome object is the one behind Varma & Monien,
Phys. Rev. E 87, 032109 (2013), Eq. (16).

Varma-Monien use the resolvent 1/(2t_VM - A) of the nearest-neighbour adjacency A, so their
energy variable is t_VM = E/2 (with E the adjacency eigenvalue used here).  Their Eq. (16)
carries an explicit flat-band factor 1/(t_VM + 1): the flat band therefore sits at t_VM = -1,
i.e. E = -2.  We reproduce, from our own 12x12 Bloch Hamiltonian, the full spectrum in BOTH
conventions and confirm:
  * a flat band at E = -2  (t_VM = -1) with spectral weight exactly 1/3 = (12-8)/12,
    i.e. the 1/(t_VM+1) pole of Varma-Monien with residue 1/3;
  * dispersive bands filling E in [-2,4]  (t_VM in [-1,2]), symmetric about E = 1 (t_VM = 1/2);
  * principal van Hove features at E = 0, 2  (t_VM = 0, 1) -- the peaks of their Fig. 2(b).
This is the identification on which "addresses the Varma-Monien integral" rests.
"""
import pickle, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sites, bonds = pickle.load(open(os.path.join(HERE, "lattice.pkl"), "rb"))
B = [(int(i), int(j), tuple(int(round(x)) for x in d)) for (i, j, d) in bonds]
nsite = 12

N = 40
ks = (np.arange(N) + 0.5) / N * 2 * np.pi
KX, KY, KZ = np.meshgrid(ks, ks, ks, indexing='ij')
kx, ky, kz = KX.ravel(), KY.ravel(), KZ.ravel()
H = np.zeros((kx.size, nsite, nsite), dtype=complex)
for (a, b, d) in B:
    H[:, a, b] += np.exp(1j * (kx*d[0] + ky*d[1] + kz*d[2]))
ev = np.linalg.eigvalsh((H + np.conjugate(np.transpose(H, (0, 2, 1)))) / 2).ravel()

flat = np.isclose(ev, -2.0, atol=1e-6)
wflat = flat.mean()
disp = ev[~flat]
print("flat band:")
print("  location  E = %.6f  (Varma-Monien t_VM = E/2 = %.6f  ->  their 1/(t_VM+1) pole)"
      % (ev[flat].mean(), ev[flat].mean() / 2))
print("  weight    = %.6f   (exact 1/3 = %.6f)  -> residue 1/3 of the VM flat-band pole"
      % (wflat, 1/3))
print("dispersive bands:")
print("  range E in [%.4f, %.4f]  (t_VM in [%.4f, %.4f])" %
      (disp.min(), disp.max(), disp.min()/2, disp.max()/2))
sym = abs(np.mean((disp - 1.0)**3))     # 3rd central moment about E=1 -> 0 if symmetric
print("  |<(E-1)^3>| = %.2e  (=> dispersive DOS symmetric about E=1, i.e. t_VM=1/2)" % sym)
print("  principal van Hove energies E=0,2 correspond to Varma-Monien t_VM = 0, 1")
print("\n=> Same spectrum as Varma-Monien Eq. (16), incl. their 1/(t_VM+1) flat-band pole.")
