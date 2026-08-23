"""
Closed-walk moment generator: rebuilds numerics/moments230.json from numerics/lattice.pkl.

Provenance only. Nothing in the certification battery runs this script; the certified
artifacts read the committed moments230.json directly. The enumeration is deterministic,
so re-running it reproduces that file exactly.
"""
import pickle, json, os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sites, bonds = pickle.load(open(os.path.join(HERE, "lattice.pkl"), "rb"))
B = [(int(i), int(j), tuple(int(round(x)) for x in d)) for (i, j, d) in bonds]
adj = defaultdict(list)
for i, j, d in B:
    adj[i].append((j, d))

NMAX = 230
st = defaultdict(int); st[(0, 0, 0, 0)] = 1; mom = [1]
for n in range(1, NMAX + 1):
    nx = defaultdict(int)
    for (b, x, y, z), c in st.items():
        for (jj, (dx, dy, dz)) in adj[b]:
            nx[(jj, x + dx, y + dy, z + dz)] += c
    st = nx
    mom.append(st.get((0, 0, 0, 0), 0))
    if n % 30 == 0:
        print("  m%d  (%d states)" % (n, len(st)), flush=True)
json.dump(mom, open(os.path.join(HERE, "moments230.json"), "w"))
print("wrote moments230.json  (m0..m%d)" % NMAX)
