import tilingdata.wakatsuki
import tilingdata.t363243432434
import sympy

from a_cyc import get_cyc
from b_p import get_p
from _1_cs import get_cs
from _3_C1 import get_C1, get_C1_undirected
from _5_C2 import get_C2, get_C2_strict_undirected
from invariants_d2 import invariants_d2
from gf_d2 import get_GF
from strict_C2 import get_C2_strict
from util_d2 import *
from cycpoly import cyccl, get_R

data = tilingdata.wakatsuki
# data = tilingdata.t363243432434
print("# V of quotient graph".ljust(30),":",data.c)

is_UD = is_undirected(data)
print("is_undirected".ljust(30),":",is_UD)

x0=0
# x0=4
cs = get_cs(data,x0,10)
print("coordination sequence".ljust(30),":",cs)

cyc = get_cyc(data)
P = get_p(data,cyc)
print("# V of P".ljust(30),":",len(P.vertices))

if is_UD: C1 = get_C1_undirected(data,x0,P)
else: C1 = get_C1(data,x0,P)
print("C1".ljust(30),":",C1,"=",C1.evalf())

C2 = get_C2(data,x0,cyc,P)
print("C2'".ljust(30),":",C2,"=",C2.evalf())

beta, denom, cpx, beta2, denom_pair= invariants_d2(data,x0,cyc,P,C1,C2)
print("beta".ljust(30),":",beta,"=",beta.evalf())
print("beta'".ljust(30),":",beta2,"=",beta2.evalf())
print("cpx".ljust(30),":",cpx)
#print("exponents of denom.".ljust(30),":",denom)
sympy.var('t')
R = get_R(denom_pair,t)
deg_R = sympy.degree(R,t)
print("degree of R".ljust(30),":",deg_R)
R_ = cyccl(R,t)
print("R'".ljust(30),":",R_)

gf = get_GF(data,x0,beta,denom)
# print("numerator".ljust(30),":",gf[0])
# print("denominator".ljust(30),":",gf[1])
print("generating function".ljust(30),":",gf[2])
temp = cyccl(gf[1],t)
denom_=1
for d in temp: denom_=denom_*(1-t**d)
numer_=sympy.polys.polytools.quo(gf[0]*denom_,gf[1])
print("numerator'".ljust(30),":",numer_)
print("denominator'".ljust(30),":",denom_)

print("Do you need the exact value of C2?(y/n)")
res = input()
if res=='y':
    C2 = get_C2_strict_undirected(data,x0,P,beta2)
    print("exact C2".ljust(30),":",C2,"=",C2.evalf())
