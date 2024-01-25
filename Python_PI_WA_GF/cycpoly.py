import sympy

def cyccl(poly,t):
    N=2*sympy.degree(poly)**2 
    fl=sympy.factor_list(poly,t)
    cl=[]
    newcl=[]
    for i in range(len(fl[1])):
        d=1
        while d < N+1 and sympy.polys.polytools.rem((1-t**d),fl[1][i][0])!=0: d+=1#
        else:cl+=[d,]*fl[1][i][1]
    while len(cl)!=0:
        n=max(cl)
        newcl.append(n)
        for i in range(1,n+1):
            if n%i==0 and i in cl: cl.remove(i)
    return(newcl)
    
def get_R(denom_pair,t):
    f = sympy.Integer(1)
    for dp in denom_pair:
        g = (1-t**dp[0]) * (1-t**dp[1])
        f = sympy.lcm(f,g)
    return f
