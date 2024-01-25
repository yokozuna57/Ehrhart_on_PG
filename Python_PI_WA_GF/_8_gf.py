import sympy

from _1_cs import get_cs,get_cs_undirected
from a_cyc import get_cyc
from b_p import get_p,mu,d_P_Gamma

sympy_t=sympy.symbols('t')

def get_GF(data,s,denom):
    den = 1
    for d in denom:
        den=den*(1-sympy_t**d)
    # print(den)
    deg = sympy.degree(den)
    cs = get_cs_undirected(data,s,deg+1)
    # cs = get_cs(data,s,deg+1)
    f = 0
    for i in range(deg+1):
        f = f + cs[i] * sympy_t**i
    f = f*den
    f = sympy.polys.polytools.rem(f,sympy_t**(deg+1))
    gc = sympy.gcd(f,den)
    f = sympy.polys.polytools.quo(f,gc)
    den = sympy.polys.polytools.quo(den,gc)
    return (f,den,f/den)
