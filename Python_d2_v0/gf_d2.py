import sympy

from _1_cs import get_cs,get_cs_undirected
from a_cyc import get_cyc
from b_p import get_p
from _3_C1 import get_C1
from _5_C2 import get_C2
from b_p import mu,d_P_Gamma

sympy_t=sympy.symbols('t')

def get_GF(data,s,beta,denom):
    den = 1
    for d in denom:
        den=den*(1-sympy_t**d)
    deg=sympy.floor(beta)+1+sympy.degree(den)
    cs = get_cs_undirected(data,s,deg)
    f = 0
    for i in range(deg+1):
        f = f + cs[i] * sympy_t**i
    f = f*den
    f = sympy.polys.polytools.rem(f,sympy_t**(deg+1))
    for d in denom:
        if sympy.polys.polytools.rem(f,(1-sympy_t**d))==0:
            den=den/(1-sympy_t**d)
            f=sympy.polys.polytools.pquo(f,(1-sympy_t**d))
    ret_numer = f
    ret_denom = den
    gc = sympy.gcd(f,den)
    f = sympy.polys.polytools.quo(f,gc)
    den = sympy.polys.polytools.quo(den,gc)
    return (ret_numer,ret_denom,f/den)
