import numpy as np
import sympy
from scipy.spatial import ConvexHull

def get_p(data,cyc):
    def mu_(x): return tuple(x[1][i]/x[0] for i in range(data.dim))
    im_mu=np.array([mu_(x) for x in cyc])
    P=ConvexHull(im_mu)
    return P

def mu(x):
    d=len(x[1])
    ret=sympy.Matrix(d,1,[0,]*d)
    for j in range(d):
        ret[j]=sympy.Rational(x[1][j],x[0])
    return ret
    
def get_vp(data,cyc,P):
    ret=[mu(cyc[i]) for i in P.vertices]
    return ret
    
def d_P_Gamma(x,cyc,P):
    for simplex in P.simplices:
        M=sympy.Matrix(len(simplex),0,[])
        for i in range(len(simplex)):
            M=M.col_insert(i,mu(cyc[simplex[i]]))
        a=M.inv()*x
        ret=[True,sympy.Integer(0)]
        for i in range(len(simplex)):
            if a[i]<0: ret[0]=False
            else: ret[1]+=a[i]
        if ret[0]: return ret[1]
