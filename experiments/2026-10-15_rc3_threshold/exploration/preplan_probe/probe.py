import sys, math
from fractions import Fraction as Q
# public export: the absolute path of the working directory, replaced by the repository root
ROOT=str(__import__('pathlib').Path(__file__).resolve().parents[4])
sys.path.insert(0, ROOT+'/experiments/2026-10-08_neumann_14rows/code')
import neumann as nm
sys.path.insert(0, ROOT+"/experiments/2026-09-18_interval_certification/code"); import certify as C
from families import build as build0, dist_to_parts
import flint, numpy as np
from flint import arb, fmpq
flint.ctx.prec=128
THR = math.sqrt(2)-1

def build_w(name,k,w):
    G = build0(name,k)
    w = Q(w)
    for e in list(G['bridges']):
        G['edges'][e]=w; G['bridges'][e]=w
    G['w']=w; G['mu']=max(abs(x) for x in G['edges'].values())
    return G

def S_exact(G):
    n=G['n']; H=np.zeros((n,n))
    for v in range(n): H[v,v]=float(G['h'][v])
    for (u,v),x in G['edges'].items(): H[u,v]=H[v,u]=float(x)
    lam,Qm=np.linalg.eigh(H)
    s,t=G['s'],G['t']
    # group
    groups=[[0]]
    for j in range(1,n):
        if lam[j]-lam[j-1]<1e-9: groups[-1].append(j)
        else: groups.append([j])
    return sum(abs(sum(Qm[t,j]*Qm[s,j] for j in g)) for g in groups)

def sig(P, d):
    s,t=P['s'],P['t']
    return sum(float(nm.upper_q(arb(d[(i,s)])))*float(nm.upper_q(arb(d[(i,t)]))) for i in range(P['k']))

def probe(name,k,w,M=300,kgrid=50):
    G=build_w(name,k,w)
    nm.build=lambda nn,kk: G
    P=nm.rig_setup(name,k)
    out,cert=nm.bounds_123(P)
    forms,info=nm.n2prime(P)
    r=dict(eta=float(P['eta'].p)/float(P['eta'].q))
    r['bridge']=sig(P,out['N1']); r['walk']=sig(P,out['N2'])
    r['chain']=sig(P,forms['N2p_path']) if 'N2p_path' in forms else float('nan')
    d={}
    for i in range(k):
        rr,_=nm.n3x_upper_reduced(P,i,M); 
        for x,v in rr.items(): d[(i,x)]=v
    r['contour']=sig(P,d)
    par=C.parameters(G); geo=C.geometry(G)
    best=min(float(C.upper_q(C.sigma1(G,par,geo,fmpq(j,kgrid)))) for j in range(1,kgrid))
    r['CT']=best
    r['S']=S_exact(G)
    return r

fams={'cliques':[0.1,3.0,3.5,4.0,4.5],'dimer':[0.05,2.5,3.0,3.3],'p4':[0.3,0.33,0.36,0.40]}
kmax={'cliques':10,'dimer':10,'p4':10}
for name,ws in fams.items():
    for w in ws:
        print(f'--- {name} w={w}')
        for k in range(2,kmax[name]+1):
            try:
                r=probe(name,k,Q(str(w)))
            except Exception as e:
                print(k,'ERR',repr(e)); continue
            flag=lambda x: '*' if x<THR else ' '
            print(f"k={k:2d} eta={r['eta']:.3f} S={r['S']:.3e}{flag(r['S'])} contour={r['contour']:.3e}{flag(r['contour'])} chain={r['chain']:.3e}{flag(r['chain'])} walk={r['walk']:.3e}{flag(r['walk'])} bridge={r['bridge']:.3e}{flag(r['bridge'])} CT={r['CT']:.3e}{flag(r['CT'])}", flush=True)
