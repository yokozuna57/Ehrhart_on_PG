import tilingdata.wakatsuki
import tilingdata.t363243432434
import tilingdata.t3122
import tilingdata.cairo
import sympy

from a_cyc import get_cyc
from b_p import get_p, get_vp
from _3_C1 import get_C1, get_C1_undirected
from _4_C2_Pinitial import get_C2_PI
from util_d2 import *
from _2_Pinitial import get_Pinitials
from cycpoly import cyccl, get_R

data = sympify(tilingdata.wakatsuki)
print("# V of quotient graph".ljust(30),":",data.c)

is_UD = is_undirected(data)
print("is_undirected".ljust(30),":",is_UD)

x0=0

cyc = get_cyc(data)
P = get_p(data,cyc)
print("# V of P".ljust(30),":",len(P.vertices))

if is_UD: C1 = get_C1_undirected(data,x0,P)
else: C1 = get_C1(data,x0,P)
print("C1".ljust(30),":",C1,"=",C1.evalf())

vp = get_vp(data,cyc,P)
is_PI = x0 in get_Pinitials(data,cyc,vp)
print("is_PI".ljust(30),":",is_PI)

if is_PI:
    C2 = get_C2_PI(data,x0,cyc,P)
    print("C2".ljust(30),":",C2,"=",C2.evalf())
