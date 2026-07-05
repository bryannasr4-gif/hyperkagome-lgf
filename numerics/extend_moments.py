import pickle,json,time
from collections import defaultdict
sites,bonds=pickle.load(open("numerics/lattice.pkl","rb"))
B=[(int(i),int(j),tuple(int(round(x)) for x in d)) for (i,j,d) in bonds]
adj=defaultdict(list)
for i,j,d in B: adj[i].append((j,d))
NMAX=230; t0=time.time()
st=defaultdict(int); st[(0,0,0,0)]=1; mom=[1]
for n in range(1,NMAX+1):
    nx=defaultdict(int)
    for (b,x,y,z),c in st.items():
        for (jj,(dx,dy,dz)) in adj[b]: nx[(jj,x+dx,y+dy,z+dz)]+=c
    st=nx; mom.append(st.get((0,0,0,0),0))
    if n%30==0: print("  m%d  (%.0fs, %d states)"%(n,time.time()-t0,len(st)),flush=True)
json.dump(mom,open("numerics/moments230.json","w"))
print("DONE m230 in %.0fs"%(time.time()-t0))
