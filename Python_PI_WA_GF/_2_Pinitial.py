from b_p import mu
def get_Pinitials(data,cyc,vp):
    j=(1<<data.c)-1
    for v in vp:
        j_=0
        for c in cyc:
            if mu(c)==v: j_=j_|c[2]
        j=j&j_
    ret=[]
    for i in range(data.c):
        if ((j>>i)&1)==1: ret.append(i)
    return ret

def get_Pinitials2(data,cyc_2,vp):
    ret=[]
    for s in range(data.c):
        f=True
        for v in vp:
            f_=False
            for c in cyc_2[s]:
                if mu(c)==v:
                    f_=True
                    break
            if not f_:
                f=False
                break
        if f: ret.append(s)
    return ret
