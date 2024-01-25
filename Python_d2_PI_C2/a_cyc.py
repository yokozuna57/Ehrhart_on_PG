def get_cyc(data):
    ret=[]
    for s in range(data.c-1,-1,-1):
        vec=[set() for i in range(1<<(s+1))]
        vec[1<<s].add((s,tuple(0 for i in range(data.dim))))
        for i in range(1<<(s+1)):
            pop_cnt=bin(i).count("1")
            for p in vec[i]:
                for ed in data.edges[p[0]]:
                    if ed[0]==s:
                        nex=tuple(p[1][i]+ed[1][i] for i in range(data.dim))
                        ret.append((pop_cnt,nex,i))
                    elif ed[0]<s and ((i>>ed[0])&1)==0:
                        nex=tuple(p[1][i]+ed[1][i] for i in range(data.dim))
                        vec[i+(1<<ed[0])].add((ed[0],nex))
    return list(set(ret))
