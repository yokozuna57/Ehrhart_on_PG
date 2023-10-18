def get_cs(data,s,N):
    b=set()
    bd={(s,)+(0,)*data.dim}
    sn=[1]
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
        sn.append(len(bd))
    return sn
    
def get_cs_undirected(data,s,N):
    pre=set()
    now={(s,)+(0,)*data.dim}
    sn=[1]
    for k in range(N):
        nex=set()
        for u in now:
            for v in data.edges[u[0]]:
                x=[v[0]]
                for i in range(data.dim):
                    x.append(u[i+1]+v[1][i])
                nex.add(tuple(x))
        nex=nex-pre-now
        sn.append(len(nex))
        pre=now
        now=nex
    return sn
