import numpy as np
import sympy
from scipy.spatial import ConvexHull

def nu(x): return sympy.ImmutableMatrix(2,1,[sympy.Rational(x[1][j],x[0]) for j in range(2)])

def get_p(data,cyc):
    points=[]
    for x in cyc: points.append(nu(x))
    P=sympy.convex_hull(*points)
    return P

'''def mu(x):
    d=len(x[1])
    ret=sympy.Matrix(d,1,[0,]*d)
    for j in range(d):
        ret[j]=sympy.Rational(x[1][j],x[0])
    return ret'''
    
def d_P_Gamma(x,P):
    for i in range(len(P.vertices)):
        j = i+1
        if j==len(P.vertices): j=0
        M=sympy.Matrix(2,2,[P.vertices[i].x,P.vertices[j].x,P.vertices[i].y,P.vertices[j].y])
        a=M.inv()*x
        ret=[True,sympy.Integer(0)]
        for i in range(2):
            if a[i]<0: ret[0]=False
            else: ret[1]+=a[i]
        if ret[0]: return ret[1]

def get_vp(data,cyc,P):
    #ret=[nu(cyc[i]) for i in P.vertices]
    ret=[]
    for i in range(len(P.vertices)):
        v=P.vertices[i]
        ret.append(nu((1,(v.x,v.y))))
    return ret
