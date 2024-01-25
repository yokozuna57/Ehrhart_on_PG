# Ehrhart_on_PG

In this repository, we will distribute related materials for our papers on growth sequences (also known as coordination sequences in crystallography) of periodic graphs. The relevant papers are as follows:

[1] T. Inoue and Y. Nakamura, Ehrhart theory on periodic graphs, [arXiv: 2305.08177v1](https://arxiv.org/abs/2305.08177)

[2] T. Inoue and Y. Nakamura, Stratified Ehrhart ring theory on periodic graphs, [arXiv: 2310.19569](https://arxiv.org/abs/2310.19569)

- A program for checking P-initial and/or well-arranged conditions of periodic graphs and calculating their growth series (i.e. generating function of growth sequences) when they are well-arranged: see [Python_PI_WA_GF](./Python_PI_WA_GF).
- A program for calculating growth series of 2D periodic graphs: see [Python_d2_v0](./Python_d2_v0).
- A program for checking P-initial condition of 2D periodic graphs and calculating C1 and C2 when they are P-initial: see [Python_d2_PI_C2](./Python_d2_PI_C2).

## input format

Here, we will explain the input format shortly.
An input should be prepared as a Python source code which includes the variants `dim`, `c`, `edges` and `pos`.
They mean the dimension, the number of vertices of the quotient graph, the edges of the quotient graph and the data of periodic realization, respectively, of the periodic graph which you want to search.
For example, the following is the input data for Wakatsuki graph ([1: Example3.9]):
```
dim=2
c=3
edges=[
    [(1,(0,0)),(1,(-1,0)),(1,(-1,-1)),(2,(0,0))],   # the targets of the edges from (0,(0,0))
    [(0,(0,0)),(0,(1,0)),(0,(1,1)),(2,(0,0))],      # the targets of the edges from (1,(0,0))
    [(0,(0,0)),(1,(0,0))]                           # the targets of the edges from (2,(0,0))
]
pos=[(0,0),(0.5,0.5),(0.5,0)]
```
(This example is [here](./Python_d2_v0/tilingdata/wakatsuki.py))

For more detail, see [2: Appendix A].
