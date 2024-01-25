import sympy
from b_p import d_P_Gamma
    
def get_C1(data,s,P):
    N=data.c-1
    b=set()
    bd={(s,)+(0,)*data.dim}
    C1=0
    for k in range(N):
        b.update(bd)
        bdk=set()
        for u in bd:
            for v in data.edges[u[0]]:
                x=[v[0]]
                for i in range(data.dim):
                    x.append(u[i+1]+v[1][i])
                bdk.add(tuple(x))
        bd=bdk.difference(b)
        for v in bd:
            C1=max(C1,d_P_Gamma(sympy.Matrix(data.dim,1,[v[i+1]+data.pos[v[0]][i]-data.pos[s][i] for i in range(data.dim)]),P)-(k+1))
    return C1

def get_C1_undirected(data,s,P):
    N=data.c-1
    pre=set()
    now={(s,)+(0,)*data.dim}
    C1=0
    for k in range(N):
        nex=set()
        for u in now:
            for v in data.edges[u[0]]:
                x=[v[0]]
                for i in range(data.dim):
                    x.append(u[i+1]+v[1][i])
                nex.add(tuple(x))
        nex=nex-pre-now
        for v in nex:
            C1=max(C1,d_P_Gamma(sympy.Matrix(data.dim,1,[v[i+1]+data.pos[v[0]][i]-data.pos[s][i] for i in range(data.dim)]),P)-(k+1))
        pre=now
        now=nex
    return C1
