import sympy

def invariants_d2(data,x0,cyc_supp,P,C1,C2):
    def nu(y): return tuple(y[1][i]/y[0] for i in range(data.dim))
    cpx=sympy.Integer(1)
    beta=0.
    beta2=0. # c.f. Theorem 4.25.
    num={}
    Len={}
    for y in cyc_supp:
        y_=sympy.ImmutableMatrix(data.dim,1,[sympy.Rational(y[1][j],y[0]) for j in range(data.dim)])
        if (y_,y[0]) not in num.keys(): num[(y_,y[0])]={y[2]}
        else: num[(y_,y[0])].add(y[2])
        if y_ not in Len.keys(): Len[y_]={y[0]}
        else: Len[y_].add(y[0])
    for x in num.keys():
        num[x]=len(num[x])
    denom=[]
    denom_pair=[]
    loop_cnt=0
    for ii in range(len(P.vertices)):
        jj = ii+1
        if jj==len(P.vertices): jj=0
        A=sympy.Matrix(2,2,[P.vertices[ii].x,P.vertices[jj].x,P.vertices[ii].y,P.vertices[jj].y])
        A=A.inv()
        V=[] # im(nu) \cap sigma
        cyc_S=[]
        for y in cyc_supp:
            y_=sympy.ImmutableMatrix(data.dim,1,[sympy.Rational(y[1][j],y[0]) for j in range(data.dim)])
            if (sympy.ones(1,data.dim)*A*y_)[0]==1:
                V.append((y_,(A*y_)[1]))
                cyc_S.append(y+((A*y_)[1],))
        V=list(set(V))
        V.sort(key=lambda x: x[1])
        m={}
        cpxs=[]
        for v in V:
            if v[1]==0 or v[1]==1:
                local_cpx=sympy.Integer(1)
                for i in range(0,len(cyc_S)):
                    if cyc_S[i][3]==v[1]: local_cpx=sympy.lcm(local_cpx,cyc_S[i][0])
                cpx=sympy.lcm(cpx,local_cpx)
                m[(v[1],v[1])]=local_cpx
                cpxs.append(local_cpx)
                continue
            local_cpx=sympy.Integer(1)
            for i in range(0,len(cyc_S)):
                for j in range(i+1,len(cyc_S)):
                    A=sympy.Matrix(2,2,[
                            sympy.Rational(cyc_S[i][1][0],cyc_S[i][0]), sympy.Rational(cyc_S[j][1][0],cyc_S[j][0]),
                            sympy.Rational(cyc_S[i][1][1],cyc_S[i][0]), sympy.Rational(cyc_S[j][1][1],cyc_S[j][0]),
                        ]
                    )
                    if A.det()==0: continue
                    A=A.inv()
                    a=A*v[0]
                    if a[0]<0 or a[1]<0: continue
                    local_cpx=sympy.lcm(local_cpx,sympy.lcm((a[0]/cyc_S[i][0]).q,(a[1]/cyc_S[j][0]).q))
            cpx=sympy.lcm(cpx,local_cpx)
            cpxs.append(local_cpx)
            for i in range(0,len(cyc_S)):
                for j in range(i+1,len(cyc_S)):
                    A=sympy.Matrix(2,2,[
                            sympy.Rational(cyc_S[i][1][0],cyc_S[i][0]), sympy.Rational(cyc_S[j][1][0],cyc_S[j][0]),
                            sympy.Rational(cyc_S[i][1][1],cyc_S[i][0]), sympy.Rational(cyc_S[j][1][1],cyc_S[j][0]),
                        ]
                    )
                    if A.det()==0: continue
                    A=A.inv()
                    a=A*v[0]
                    if a[0]<0 or a[1]<0: continue
                    upd=[
                        ((v[1],cyc_S[i][3]),a[0]*local_cpx),
                        ((v[1],cyc_S[j][3]),a[1]*local_cpx),
                    ]
                    for temp in upd:
                        if temp[0] not in m.keys(): m[temp[0]]=temp[1]
                        else: m[temp[0]]=sympy.Max(m[temp[0]],temp[1])
        for x in cpxs:
            if x not in denom: denom.append(x)
        for i in range(len(cpxs)-1):
            if cpxs[i]==cpxs[i+1]:
                if denom.count(x)==1: denom.append(x)
        for i in range(len(cpxs)-1):
            denom_pair.append((cpxs[i],cpxs[i+1]))
        for i in range(len(V)-1):
            # delta = {V[i],V[i+1]}
            beta_delta=sympy.Integer(0)
            beta2_delta=sympy.Integer(0)
            # v=V[i]
            v=V[i]
            beta_v=sympy.Integer(0)
            # v=V[i], F={V[i],V[i+1],...}
            if i>0:
                F=[V[j][0] for j in range(i,len(V))]
                s=0
                for j in range(i,len(V)):
                    u=V[j]
                    for d in Len[u[0]]:
                        s+=m[(v[1],u[1])]+d*(num[(u[0],d)]-1)
                h=sympy.Rational(1,1000)
                a=sympy.Rational(1,1000)
                A=V[i-1][0].col_insert(1,V[len(V)-1][0])
                A=A.inv()
                for y in cyc_supp:
                    y_=sympy.ImmutableMatrix(data.dim,1,[sympy.Rational(y[1][j],y[0]) for j in range(data.dim)])
                    if y_ not in F:
                        aa=A*y_
                        if aa[1]>0: h=sympy.Max(h,aa[1]/(1-aa[0]))
                        # note that u = aa[0]V[i-1] + aa[1]V[len-1] = aa[0]V[i-1] + (1-aa[0]) aa[1]/(1-aa[0]) V[len-1]
                A=V[i-1][0].col_insert(1,V[len(V)-1][0]*h)
                A=A.inv()
                a=1/(sympy.ones(1,data.dim)*A*v[0])[0]
                beta_v=sympy.Max(beta_v,a/(1-a)*(C1/h+C2+(1-h)*(s+data.c-1)/h)-cpxs[i])
                # note that alpha <= alpha' - cpx because s >= m(v,v) = cpx(sigma,v) and a >= h
            # v=V[i], F={V[0],V[1],...,V[i]}
            F=[V[j][0] for j in range(0,i+1)]
            s=0
            for j in range(0,i+1):
                u=V[j]
                for d in Len[u[0]]:
                    s+=m[(v[1],u[1])]+d*(num[(u[0],d)]-1)
            h=sympy.Rational(1,1000)
            a=sympy.Rational(1,1000)
            A=V[i][0].col_insert(1,V[i+1][0])
            A=A.inv()
            for y in cyc_supp:
                y_=sympy.ImmutableMatrix(data.dim,1,[sympy.Rational(y[1][j],y[0]) for j in range(data.dim)])
                if y_ not in F:
                    aa=A*y_
                    if aa[0]>0: a=sympy.Max(a,aa[0]/(1-aa[1]))
            A=(V[i][0]*a).col_insert(1,V[i+1][0])
            A=A.inv()
            h=1/(sympy.ones(1,data.dim)*A*V[0][0])[0]
            beta_v=sympy.Max(beta_v,a/(1-a)*(C1/h+C2+(1-h)*(s+data.c-1)/h)-cpxs[i])
            beta_delta+=beta_v
            beta2_delta+=beta_v+cpxs[i]
            # v=V[i+1]
            v=V[i+1]
            beta_v=sympy.Integer(0)
            # v=V[i+1], F={V[0],V[1],...,V[i+1]}
            if i+1<len(V)-1:
                F=[V[j][0] for j in range(0,i+2)]
                s=0
                for j in range(0,i+2):
                    u=V[j]
                    for d in Len[u[0]]:
                        s+=m[(v[1],u[1])]+d*(num[(u[0],d)]-1)
                h=sympy.Rational(1,1000)
                a=sympy.Rational(1,1000)
                A=V[i+2][0].col_insert(1,V[0][0])
                A=A.inv()
                for y in cyc_supp:
                    y_=sympy.ImmutableMatrix(data.dim,1,[sympy.Rational(y[1][j],y[0]) for j in range(data.dim)])
                    if y_ not in F:
                        aa=A*y_
                        if aa[1]>0: h=sympy.Max(h,aa[1]/(1-aa[0]))
                A=V[i+2][0].col_insert(1,V[0][0]*h)
                A=A.inv()
                a=1/(sympy.ones(1,data.dim)*A*v[0])[0]
                beta_v=sympy.Max(beta_v,a/(1-a)*(C1/h+C2+(1-h)*(s+data.c-1)/h)-cpxs[i+1])
            # v=V[i+1], F={V[i+1],V[i+2],...}
            F=[V[j][0] for j in range(i+1,len(V))]
            s=0
            for j in range(i+1,len(V)):
                u=V[j]
                for d in Len[u[0]]:
                    s+=m[(v[1],u[1])]+d*(num[(u[0],d)]-1)
            h=sympy.Rational(1,1000)
            a=sympy.Rational(1,1000)
            A=V[i+1][0].col_insert(1,V[i][0])
            A=A.inv()
            for y in cyc_supp:
                y_=sympy.ImmutableMatrix(data.dim,1,[sympy.Rational(y[1][j],y[0]) for j in range(data.dim)])
                if y_ not in F:
                    aa=A*y_
                    if aa[0]>0: a=sympy.Max(a,aa[0]/(1-aa[1]))
            A=(V[i+1][0]*a).col_insert(1,V[i][0])
            A=A.inv()
            h=1/(sympy.ones(1,data.dim)*A*V[len(V)-1][0])[0]
            beta_v=sympy.Max(beta_v,a/(1-a)*(C1/h+C2+(1-h)*(s+data.c-1)/h)-cpxs[i+1])
            beta_delta+=beta_v
            beta2_delta+=beta_v+cpxs[i+1]
            beta=sympy.Max(beta,beta_delta)
            beta2=sympy.Max(beta2,beta2_delta)
    beta=beta+C2
    beta2=beta2+C2
    return (beta,denom,cpx,beta2,denom_pair)
