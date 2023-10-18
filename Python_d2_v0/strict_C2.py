import networkx as nx
import sympy
import math
from b_p import d_P_Gamma

def get_C2_strict(data,s,cyc,P,beta2): #networkxを使わないように書き直したいかも
    MAX_d=math.floor(beta2)
    MIN_x=tuple(min([e[1][j] for i in range(data.c) for e in data.edges[i]])*MAX_d for j in range(data.dim))
    MAX_x=tuple(max([e[1][j] for i in range(data.c) for e in data.edges[i]])*MAX_d for j in range(data.dim))
    G=nx.DiGraph()
    x=list(MIN_x)
    while True:
        G.add_edges_from([((i,tuple(x)),(e[0],tuple(x[j]+e[1][j] for j in range(data.dim))))
        for i in range(data.c) for e in data.edges[i]])
        if x==list(MAX_x):
            break
        for i in range(data.dim-1,-1,-1):
            if x[i]!=MAX_x[i]:
                x[i]+=1
                break
            else:
                x[i]=MIN_x[i]
    dist=nx.shortest_path_length(G,source=(s,(0,)*data.dim))
    C1=0
    for v in dist:
        if dist[v]<=data.c-1:
            C1=max(C1,d_P_Gamma(sympy.Matrix(data.dim,1,[v[1][i]+data.pos[v[0]][i]-data.pos[s][i] for i in range(data.dim)]),cyc,P)-dist[v])
    return C1
