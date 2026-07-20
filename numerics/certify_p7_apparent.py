"""Certify that the degree-seven locus p7 consists of APPARENT singularities of M.

The exponents at each root of p7 are the integers {0,1,3} (computed exactly by reduction mod p7
in certify_factor.py). Integer exponents are NECESSARY but NOT SUFFICIENT for apparentness: a
point with integer exponents can still carry a logarithm. Apparentness requires in addition that
all three local solutions are log-free, i.e. that the Frobenius obstruction matrix has full
nullity 3.

The roots of p7 are algebraic of degree 7, so the check is done over GF(p) for primes p at which
p7 has roots; the seven roots are Galois-conjugate, so a single root already decides the question,
but every available root at every tested prime is checked.

Validated first on two operators with known structure:
    t^3 D^3  -> 3 log-free solutions (apparent, no log)
    theta^3  -> 1 log-free solution  (genuine log)

Run:  python certify_p7_apparent.py      (exit 0 on success)
"""
import json
import os
import sys
from collections import defaultdict

import sympy as sp

t = sp.symbols('t')
s = sp.symbols('s')

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, 'M_coeffs.json')) as fh:
    M_COEFFS = [[int(a) for a in row] for row in json.load(fh)['coeffs']]

# p7 in ASCENDING order of powers of t
P7_ASC = [-64, 1800, -17580, 93657, -290956, 455798, -369600, 101025]


def gf_rank(rows, p):
    """Rank of an integer matrix over GF(p) by Gaussian elimination."""
    m = [r[:] for r in rows]
    nrows, ncols = len(m), len(m[0])
    rank = 0
    for c in range(ncols):
        piv = next((r for r in range(rank, nrows) if m[r][c] % p), None)
        if piv is None:
            continue
        m[rank], m[piv] = m[piv], m[rank]
        inv = pow(m[rank][c], p - 2, p)
        m[rank] = [(v * inv) % p for v in m[rank]]
        for r in range(nrows):
            if r != rank and m[r][c] % p:
                f = m[r][c]
                m[r] = [(m[r][j] - f * m[rank][j]) % p for j in range(ncols)]
        rank += 1
    return rank


def count_log_free_modp(coeffs_asc, pt, p, K=10):
    """(#log-free local solutions at t=pt over GF(p), small integer exponents there)."""
    n = len(coeffs_asc) - 1
    # shift each coefficient polynomial to the local variable x = t - pt, over GF(p)
    cij = []
    for i in range(n + 1):
        co = [c % p for c in coeffs_asc[i]]
        shifted = [0] * len(co)
        for j, c in enumerate(co):                      # c * (pt + x)^j
            if c == 0:
                continue
            for k in range(j + 1):
                shifted[k] = (shifted[k] + c * int(sp.binomial(j, k)) * pow(pt, j - k, p)) % p
        cij.append({j: v for j, v in enumerate(shifted) if v})

    # q_d(s) = sum_i c_{i,j} * s(s-1)...(s-i+1) collected by d = j - i
    qd = defaultdict(lambda: [0])

    def addpoly(a, b):
        out = [0] * max(len(a), len(b))
        for idx, v in enumerate(a):
            out[idx] = (out[idx] + v) % p
        for idx, v in enumerate(b):
            out[idx] = (out[idx] + v) % p
        return out

    for i in range(n + 1):
        fall = sp.Poly(sp.prod([s - a for a in range(i)]), s) if i else sp.Poly(1, s)
        fc = [int(c) % p for c in reversed(fall.all_coeffs())]
        for j, c in cij[i].items():
            qd[j - i] = addpoly(qd[j - i], [(c * v) % p for v in fc])

    def ev(poly, val):
        acc = 0
        for c in reversed(poly):
            acc = (acc * val + c) % p
        return acc

    dmin = min(d for d in qd if any(qd[d]))
    exps = [e for e in range(0, 25) if ev(qd[dmin], e) == 0]   # exponents here are small integers
    emin = 0
    rows = []
    for r in range(K):
        rows.append([ev(qd[dmin + r - k], emin + k) if (k <= r and (dmin + r - k) in qd) else 0
                     for k in range(K)])
    return K - gf_rank(rows, p), exps


def main():
    print('Validating the mod-p log-free counter on operators with known structure:')
    v_nolog, _ = count_log_free_modp([[0], [0], [0], [0, 0, 0, 1]], 0, 10**9 + 7)   # t^3 D^3
    v_log, _ = count_log_free_modp([[0], [0, 1], [0, 0, 3], [0, 0, 0, 1]], 0, 10**9 + 7)  # theta^3
    print(f'  t^3 D^3 : log-free = {v_nolog} of 3  (expect 3, apparent)')
    print(f'  theta^3 : log-free = {v_log} of 3  (expect 1, genuine log)')
    assert (v_nolog, v_log) == (3, 1), 'mod-p log-free counter FAILED validation'

    p7poly = sp.Poly(list(reversed(P7_ASC)), t)
    print('\nTesting every root of p7 available at each of several primes:')
    tested_primes, tested_roots, ok = 0, 0, True
    for p in sp.primerange(10007, 60000):
        if any(c % p == 0 for c in (101025, 64)):
            continue
        roots = sorted(int(r) % p for r in sp.Poly(p7poly, t, modulus=p).ground_roots())
        if not roots:
            continue
        print(f'  prime p = {p}: p7 has {len(roots)} root(s) mod p')
        for r in roots:
            nlf, exps = count_log_free_modp(M_COEFFS, r, p)
            good = (nlf == 3)
            ok &= good
            print(f'     root {r:>8}: log-free = {nlf} of 3, exponents {exps}'
                  f'  -> {"APPARENT" if good else "NOT APPARENT (genuine log)"}')
            tested_roots += 1
        tested_primes += 1
        if tested_primes == 5:
            break

    print()
    if ok and tested_roots:
        print(f'RESULT: PASS. All {tested_roots} tested roots of p7 (over {tested_primes} primes) have')
        print('  3 log-free local solutions of 3 => p7 is a locus of APPARENT singularities of M,')
        print('  with exponents {0,1,3}. (Integer exponents alone would not establish this.)')
        return 0
    print('RESULT: FAIL -- a root of p7 carries a logarithm; it is not apparent.')
    return 1


if __name__ == '__main__':
    sys.exit(main())
