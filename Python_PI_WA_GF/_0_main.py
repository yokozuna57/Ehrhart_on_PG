import tilingdata.wakatsuki
import tilingdata.cairo
import tilingdata.snub632
import tilingdata.t346
import tilingdata.t3464
import tilingdata.t482
import tilingdata.t3122
import tilingdata.t32434
import tilingdata.t363243432434
import tilingdata.t363243432434bai
import tilingdata.SACADA1
import tilingdata.SACADA60
import SACADApy.SACADA1p
import SACADApy.SACADA8p
import SACADApy.SACADA10p
import SACADApy.SACADA11p
import SACADApy.SACADA12p
import SACADApy.SACADA13p
import SACADApy.SACADA20p
import SACADApy.SACADA21p
import SACADApy.SACADA29p
import SACADApy.SACADA30p
import SACADApy.SACADA31p
import SACADApy.SACADA33p
import SACADApy.SACADA35p
import SACADApy.SACADA37p
import SACADApy.SACADA39p
import SACADApy.SACADA51p
import SACADApy.SACADA52p
import SACADApy.SACADA56p
import SACADApy.SACADA57p
import SACADApy.SACADA58p
import SACADApy.SACADA59p
import SACADApy.SACADA60p
import SACADApy.SACADA65p
import SACADApy.SACADA66p
import SACADApy.SACADA67p
import SACADApy.SACADA68p
import SACADApy.SACADA69p
import SACADApy.SACADA70p
import SACADApy.SACADA71p
import SACADApy.SACADA73p
import SACADApy.SACADA74p
import SACADApy.SACADA75p
import SACADApy.SACADA76p
import SACADApy.SACADA77p
import SACADApy.SACADA78p
import SACADApy.SACADA79p
import SACADApy.SACADA80p
import SACADApy.SACADA81p
import SACADApy.SACADA82p
import SACADApy.SACADA83p
import SACADApy.SACADA84p
import SACADApy.SACADA85p
import SACADApy.SACADA86p
import SACADApy.SACADA87p
import SACADApy.SACADA88p
import SACADApy.SACADA89p
import SACADApy.SACADA613p
import SACADApy.SACADA628p
import SACADApy.SACADA629p
import sympy
import sys

from _1_cs import get_cs,get_cs_undirected
from a_cyc import get_cyc,get_cyc_,get_cyc2
from b_p import get_p,get_vp
from _2_Pinitial import get_Pinitials, get_Pinitials2
from _6_WellArranged import is_WellArranged
from _8_gf import get_GF
from cycpoly import cyccl, get_R

#args = sys.argv
#exec(f"data=SACADApy.SACADA{args[1]}p")
data=SACADApy.SACADA1p
#data=tilingdata.cairo
#data=tilingdata.t32434
#data=tilingdata.t482
#data=tilingdata.t3464
#data=tilingdata.t3122
#data=tilingdata.t346
#data=tilingdata.snub632

print("# V of quotient graph".ljust(30),":",data.c)

cyc=get_cyc(data)
P=get_p(data,cyc)
print("# V of P".ljust(30),":",len(P.vertices))

vp = get_vp(data,cyc,P)
#Pinitials = get_Pinitials(data,cyc,vp)
#print("P-initial vertices".ljust(30),":",Pinitials)

cyc_2 = get_cyc2(data)
Pinitials = get_Pinitials2(data,cyc_2,vp)
print("P-initial vertices".ljust(30),":",Pinitials)

for s in Pinitials:
    print("s:",s)
    is_WA, denom = is_WellArranged(data,s,cyc,P,cyc_2)
    if is_WA:
        print("is_WA".ljust(30),":",is_WA)
        sympy.var('t')
        gf = get_GF(data,s,denom)
        print("generating function".ljust(30),":",gf[2])
        temp = cyccl(gf[1],t)
        denom_=1
        for d in temp: denom_=denom_*(1-t**d)
        numer_=sympy.polys.polytools.quo(gf[0]*denom_,gf[1])
        print("numerator'".ljust(30),":",numer_)
        print("denominator'".ljust(30),":",denom_)
    else:
        print("is_WA".ljust(30),":","unknown")
        print("GF".ljust(30),":","-")
