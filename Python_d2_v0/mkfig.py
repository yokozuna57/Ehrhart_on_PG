import numpy as np
import matplotlib.pyplot as plt

def graphview(tdata):
    N=5
    x=[((np.array([a0,a1])+np.array(tdata.pos[k]))*np.matrix(tdata.lattice_const))[0,0] for a0 in range(-N,N+1) for a1 in range(-N,N+1) for k in range(tdata.c)] 
    y=[((np.array([a0,a1])+np.array(tdata.pos[k]))*np.matrix(tdata.lattice_const))[0,1] for a0 in range(-N,N+1) for a1 in range(-N,N+1) for k in range(tdata.c)]

    #plt.grid(True)
    plt.scatter(x, y, s=8, color="black") #散布図の描画

    for a0 in range(-N,N+1):
        for a1 in range(-N,N+1):
            for k in range(tdata.c):
                for u in tdata.edges[k]:
                    if u[0]>k:
                        x0=((np.array([a0,a1])+np.array(tdata.pos[k]))*np.matrix(tdata.lattice_const))[0,0]
                        y0=((np.array([a0,a1])+np.array(tdata.pos[k]))*np.matrix(tdata.lattice_const))[0,1]
                        x1=((np.array([a0,a1])+np.array(u[1])+np.array(tdata.pos[u[0]]))*np.matrix(tdata.lattice_const))[0,0]
                        y1=((np.array([a0,a1])+np.array(u[1])+np.array(tdata.pos[u[0]]))*np.matrix(tdata.lattice_const))[0,1]
                        x=(x0,x1)
                        y=(y0,y1)
                        #plt.plot(x,y,marker='.',markersize=5,color='black',lw=1)
                        plt.plot(x,y,color='black',lw=0.5)

    
    v00=(np.array([0,0])*np.matrix(tdata.lattice_const))#単位格子の描画
    v01=(np.array([0,1])*np.matrix(tdata.lattice_const))
    v10=(np.array([1,0])*np.matrix(tdata.lattice_const))
    v11=(np.array([1,1])*np.matrix(tdata.lattice_const))
    x=(v00[0,0],v10[0,0],v11[0,0],v01[0,0],v00[0,0])
    y=(v00[0,1],v10[0,1],v11[0,1],v01[0,1],v00[0,1])
    plt.plot(x,y,color='red',lw=1.5)

    for i in range(tdata.c):
        xi=(np.array([0,0])+np.array(tdata.pos[i]))*np.matrix(tdata.lattice_const)
        x=(xi[0,0],)
        y=(xi[0,1],)
        moji='x'+str(i)
        plt.scatter(x,y,s=20,color='blue')
        plt.annotate(moji, (xi[0,0]+0.1, xi[0,1]-0.1), fontsize=8, alpha=0.9)
            
    plt.xlim(-2,8)
    plt.ylim(-8,4)
    
    plt.show()


def uvecplot2d(tdata):
    P=get_p(tdata)
    ptnum=len(P.points)

    x=[P.points[i][0] for i in range(ptnum)] 
    y=[P.points[i][1] for i in range(ptnum)]
 
    plt.grid(True)
    plt.scatter(x, y, s=25, color="blue") #散布図の描画

    for i in range(ptnum): #頂点にラベルを付ける
        plt.annotate(i, (x[i], y[i]), fontsize=8)

    plt.show()


import tilingdata.cairo
import tilingdata.t32434
import tilingdata.t482
import tilingdata.t346
import tilingdata.snub632
import tilingdata.t3464
import tilingdata.t3122
#import tilingdata.t363243432434
#import tilingdata.t363243432434bai
#import tilingdata.a250125
#import tilingdata.a307270

tdata=tilingdata.cairo

graphview(tdata)




