dim=2
c=6
edges=[
    [(1,(0,0)),(2,(-1,0)),(3,(0,-1))], #0
    [(0,(0,0)),(2,(0,0)),(3,(0,0))], #1
    [(1,(0,0)),(4,(0,0)),(5,(0,-1)),(0,(1,0))], #2
    [(1,(0,0)),(5,(0,0)),(4,(-1,0)),(0,(0,1))], #3
    [(2,(0,0)),(5,(0,0)),(3,(1,0))], #4
    [(3,(0,0)),(4,(0,0)),(2,(0,1))], #5
]
pos=[
    (0,0),#0
    (1/6,1/6),#1
    (7/12,1/12),#2
    (1/12,7/12),#3
    (2/3,1/2),#4
    (1/2,2/3),#5
]
lattice_const=[
    (2,-2),
    (2,2),
]
