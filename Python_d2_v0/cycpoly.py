import sympy

# sympy.var('t')

'''f = -t**18 - 2*t**17 - 3*t**16 - 3*t**15 - t**14 + t**13 + 3*t**12 + 3*t**11 + 2*t**10 + 2*t**9 + 2*t**8 + 3*t**7 + 3*t**6 + t**5 - t**4 - 3*t**3 - 3*t**2 - 2*t - 1

print(sympy.factor_list(f,t))'''

def cyccl(poly,t):
    N=2*sympy.degree(poly)**2 #polyが絶対値1以外の根を持つ場合にループがおきるのを防止. 2(nのオイラー数)^2 >= n。
    fl=sympy.factor_list(poly,t)
    cl=[]
    newcl=[]
    for i in range(len(fl[1])):
        d=1
        while d < N+1 and sympy.polys.polytools.rem((1-t**d),fl[1][i][0])!=0: d+=1#
        else:cl+=[d,]*fl[1][i][1]
    #print(cl)#
    while len(cl)!=0:
        n=max(cl)
        newcl.append(n)
        for i in range(1,n+1):
            if n%i==0 and i in cl: cl.remove(i)
    #print(newcl)#
    return(newcl)
    
def get_R(denom_pair,t):
    f = sympy.Integer(1)
    for dp in denom_pair:
        g = (1-t**dp[0]) * (1-t**dp[1])
        f = sympy.lcm(f,g)
    return f

# cyccl(f)
