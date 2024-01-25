def get_cyc_(data):
    ret=[]
    for s in range(data.c-1,-1,-1):
        vec=[set() for i in range(1<<(s+1))]
        vec[1<<s].add((s,tuple(0 for i in range(data.dim))))
        for i in range(1<<(s+1)):
            pop_cnt=bin(i).count("1")
            for p in vec[i]:
                for ed in data.edges[p[0]]:
                    if ed[0]<=s:
                        nex=tuple(p[1][i]+ed[1][i] for i in range(data.dim))
                        if ed[0]==s: ret.append((pop_cnt,nex,i))
                        elif ((i>>ed[0])&1)==0:
                            vec[i+(1<<ed[0])].add((ed[0],nex))
    return list(set(ret))

def get_cyc(data):
    ret=[]
    for s in range(data.c):
        b=set()
        bd={(s,)+(0,)*data.dim}
        sn=[1]
        for k in range(data.c):
            b.update(bd)
            bdk=set()
            for u in bd:
                for v in data.edges[u[0]]:
                    x=[v[0]]
                    for i in range(data.dim):
                        x.append(u[i+1]+v[1][i])
                    bdk.add(tuple(x))
                    if x[0]==s:
                        ret.append((k+1,tuple(x[1:])))
            bd=bdk.difference(b)
            sn.append(len(bd))
    return list(set(ret))

def get_cyc2(data):
    ret=[]
    for s in range(data.c):
        ret.append([])
        b=set()
        bd={(s,)+(0,)*data.dim}
        sn=[1]
        for k in range(data.c):
            b.update(bd)
            bdk=set()
            for u in bd:
                for v in data.edges[u[0]]:
                    x=[v[0]]
                    for i in range(data.dim):
                        x.append(u[i+1]+v[1][i])
                    bdk.add(tuple(x))
                    if x[0]==s:
                        ret[s].append((k+1,tuple(x[1:])))
            bd=bdk.difference(b)
            sn.append(len(bd))
    return ret
