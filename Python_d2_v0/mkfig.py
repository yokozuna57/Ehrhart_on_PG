import numpy as np
import matplotlib.pyplot as plt

def graphview(tdata):
    N=5
    x=[((np.array([a0,a1])+np.array(tdata.pos[k]))*np.matrix(tdata.lattice_const))[0,0] for a0 in range(-N,N+1) for a1 in range(-N,N+1) for k in range(tdata.c)] 
    y=[((np.array([a0,a1])+np.array(tdata.pos[k]))*np.matrix(tdata.lattice_const))[0,1] for a0 in range(-N,N+1) for a1 in range(-N,N+1) for k in range(tdata.c)]

    #plt.grid(True)
    plt.scatter(x, y, s=8, color="black")

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

    
    v00=(np.array([0,0])*np.matrix(tdata.lattice_const))#lattice
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


import tilingdata.cairo
import tilingdata.t32434

tdata=tilingdata.cairo

graphview(tdata)




