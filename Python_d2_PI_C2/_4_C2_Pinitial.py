import sympy
#from b_p import mu,nu,d_P_Gamma
from b_p import nu,d_P_Gamma

def get_C2_PI(data,s,cyc,P,r=0):
    list_y1=[]
    list_y2=[]
    cyc_P=[]
    for c in cyc:
        if ((c[2]>>s)&1)==1: cyc_P.append((c[0],c[1]))
    cyc_P=list(set(cyc_P))
    for ii in range(len(P.vertices)):
        jj = ii+1
        if jj==len(P.vertices): jj=0
        ds=[]
        min_x=[0]*data.dim
        max_x=[0]*data.dim
        d=data.c
        for c in cyc_P:
            if nu(c)==sympy.ImmutableMatrix([P.vertices[ii].x,P.vertices[ii].y]): d=sympy.Min(d,c[0])
        ds.append(d)
        if P.vertices[ii].x<0: min_x[0]+=P.vertices[ii].x*d
        else: max_x[0]+=P.vertices[ii].x*d
        if P.vertices[ii].y<0: min_x[1]+=P.vertices[ii].y*d
        else: max_x[1]+=P.vertices[ii].y*d
        d=data.c
        for c in cyc_P:
            if nu(c)==sympy.ImmutableMatrix([P.vertices[jj].x,P.vertices[jj].y]): d=sympy.Min(d,c[0])
        ds.append(d)
        if P.vertices[jj].x<0: min_x[0]+=P.vertices[jj].x*d
        else: max_x[0]+=P.vertices[jj].x*d
        if P.vertices[jj].y<0: min_x[1]+=P.vertices[jj].y*d
        else: max_x[1]+=P.vertices[jj].y*d
        min_x=tuple(sympy.floor(min_x[i]) for i in range(data.dim))
        max_x=tuple(sympy.floor(max_x[i]) for i in range(data.dim))
        M=sympy.Matrix(2,2,[P.vertices[ii].x,P.vertices[jj].x,P.vertices[ii].y,P.vertices[jj].y])
        M=M.inv()
        x=list(min_x)
        while True:
            for i in range(data.c):
                pos=[x[j]+data.pos[i][j]-data.pos[s][j] for j in range(data.dim)]
                a=M*sympy.Matrix(data.dim,1,pos)
                ret=[True,True,0]
                for j in range(data.dim):
                    if a[j]<0 or a[j]>=r*ds[j]: ret[0]=False
                    if a[j]<0 or a[j]>=(r+1)*ds[j]: ret[1]=False
                    #ret[2]+=a[j]*ds[j]
                    ret[2]+=a[j]
                if ret[1]:
                    if ret[0]: list_y1.append(((i,tuple(x)),ret[2]))
                    else: list_y2.append(((i,tuple(x)),ret[2]))
            if x==list(max_x):
                break
            for i in range(data.dim-1,-1,-1):
                if x[i]!=max_x[i]:
                    x[i]+=1
                    break
                else:
                    x[i]=min_x[i]
    #list_y1=list(set(list_y1))
    #list_y2=list(set(list_y2))
    #print(list_y1,list_y2)
    list_y1=list(set(list_y2))
    list_y2=[]
    #print(list_y1,list_y2)
    que=[(((s,(tuple(0 for i in range(data.dim)))),(1<<s)),0)]
    q_nex=0
    used=set(que)
    dict_y1=dict()
    dict_y2=dict()
    for y in list_y1:
        if y[0]!=que[0][0][0]:
            dict_y1[y[0]]=y[1]
    for y in list_y2:
        dict_y2[y[0]]=y[1]
    C2=0
    while len(dict_y1)>0 or len(dict_y2)>0:
        q=que[q_nex]
        q_nex+=1
        for ed in data.edges[q[0][0][0]]:
            temp=q[0][1]
            if ((temp>>ed[0])&1)==0: temp+=1<<ed[0]
            nex=((ed[0],tuple(q[0][0][1][i]+ed[1][i] for i in range(data.dim))),temp)
            if nex in used: continue
            else :
                que.append((nex,q[1]+1))
                used.add(nex)
                if nex[0] in dict_y1:
                    C2=max(C2,q[1]+1-dict_y1[nex[0]])
                    del dict_y1[nex[0]]
                if nex[1]==((1<<data.c)-1) and nex[0] in dict_y2:
                    C2=max(C2,q[1]+1-dict_y2[nex[0]])
                    del dict_y2[nex[0]]
    return C2

def get_C2_strict_undirected(data,s,P,beta2):
    N=sympy.floor(beta2)
    pre=set()
    now={(s,)+(0,)*data.dim}
    C2=0
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
            C2=max(C2,(k+1)-d_P_Gamma(sympy.Matrix(data.dim,1,[v[i+1]+data.pos[v[0]][i]-data.pos[s][i] for i in range(data.dim)]),P))
        pre=now
        now=nex
    return C2
