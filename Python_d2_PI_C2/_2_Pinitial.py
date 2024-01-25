from b_p import nu
def get_Pinitials(data,cyc,vp):
    j=(1<<data.c)-1
    for v in vp:
        j_=0
        for c in cyc:
            if nu(c)==v: j_=j_|c[2]
        j=j&j_
    ret=[]
    for i in range(data.c):
        if ((j>>i)&1)==1: ret.append(i)
    return ret
