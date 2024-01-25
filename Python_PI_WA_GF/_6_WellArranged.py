import sympy
import math
from b_p import mu,d_P_Gamma

def f(x):
    ret = []
    for i in range(x.shape[0]):
        ret.append(x[i,0])
    return tuple(ret)

def is_WellArranged(data,s,cyc,P,cyc_2):
    dict_c={}
    for c in cyc_2[s]:
        mu_c = f(mu(c))
        if mu_c in dict_c:
            if c[0]<dict_c[mu_c][0]: dict_c[mu_c] = c
        else: dict_c[mu_c]=c
    list_y=[]
    list_pair=[]
    denom_dict={}
    for simplex in P.simplices:
        ds=[]
        M=sympy.Matrix(len(simplex),0,[])
        min_x=[0]*data.dim
        max_x=[0]*data.dim
        for i in range(len(simplex)):
            c=cyc[simplex[i]]
            c = dict_c[f(mu(c))]
            ds.append(c[0])
            M=M.col_insert(i,sympy.Matrix(len(simplex),1,c[1]))
            for j in range(data.dim):
                if c[1][j]<0: min_x[j]+=c[1][j]
                else: max_x[j]+=c[1][j]
        M=M.inv()
        for d in ds:
            cnt=0
            for d_ in ds:
                if d==d_: cnt = cnt+1
            if d in denom_dict:
                denom_dict[d]=max(denom_dict[d],cnt)
            else: denom_dict[d]=cnt
        min_x=tuple(math.floor(min_x[i]-1.) for i in range(data.dim))
        max_x=tuple(math.floor(max_x[i]) for i in range(data.dim))
        x=list(min_x)
        while True:
            for i in range(data.c):
                pos=[x[j]+data.pos[i][j]-data.pos[s][j] for j in range(data.dim)]
                a=M*sympy.Matrix(data.dim,1,pos)
                ret=True
                diag=[0]*data.dim
                dist=0
                for j in range(data.dim):
                    if a[j]<0 or a[j]>=1: ret=False
                    elif a[j]>0:
                        for k in range(data.dim): diag[k]+=dict_c[f(mu(cyc[simplex[j]]))][1][k]
                        dist = dist+ds[j]
                if ret:
                    list_y.append((i,tuple(x)))
                    list_y.append((i,tuple(x[j]-diag[j] for j in range(data.dim))))
                    list_y.append((s,tuple(diag)))
                    list_pair.append(((i,tuple(x)),(i,tuple(x[j]-diag[j] for j in range(data.dim))),(s,tuple(diag)),dist))
            if x==list(max_x):
                break
            for i in range(data.dim-1,-1,-1):
                if x[i]!=max_x[i]:
                    x[i]+=1
                    break
                else:
                    x[i]=min_x[i]
    list_y=set(list_y)
    que=[((s,(tuple(0 for i in range(data.dim)))),0)]
    q_nex=0
    dist=dict()
    dist[que[0][0]]=0
    if que[0][0] in list_y: list_y.discard(que[0][0])
    while len(list_y)>0:
        q=que[q_nex]
        q_nex+=1
        for ed in data.edges[q[0][0]]:
            nex=(ed[0],tuple(q[0][1][i]+ed[1][i] for i in range(data.dim)))
            if nex in dist: continue
            else:
                que.append((nex,q[1]+1))
                dist[nex]=q[1]+1
                if nex in list_y:
                    list_y.discard(nex)
    for p in list_pair:
        if dist[p[0]]+dist[p[1]]!=p[3]: return (False,[])
    denom_list=[]
    for d in denom_dict:
        for i in range(denom_dict[d]):
            denom_list.append(d)
    return (True,denom_list)

