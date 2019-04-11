General Topological Overlap Measure
===================================

This package takes an undirected unweighted scipy.sparse adjacency
matrix as an input and computes the GTOM(m) method, using m+1-step
neighbors (1). It's highly efficient and can be used for parallel
computation by calling the function for only few nodes at a time.

(1) `*Gene network interconnectedness and the generalized topological
overlap
measure* <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.93.329&rep=rep1&type=pdf>`__,
A. M. Yip and S. Horvath, BMC Bioinformatics 2007 8:22

Install
-------

.. code:: bash

    $ pip install gtom

Example
-------

Recreating Panels a and c of Figure 3 of the technical report [1]:

.. code:: python

    #!python

    import matplotlib.pyplot as pl
    import networkx as nx
    from gtom import gtom
    import scipy.sparse as sprs
    import numpy as np

    edges = [(0,1),(1,2),(0,3),(0,4),(0,5),(0,7),
             (1,3),(1,4),(1,6),(1,8),(1,9),(1,10),
             (5,6),(7,8)]
    G = nx.Graph()
    G.add_edges_from(edges)

    pos = nx.spring_layout(G)
    labels = { n:str(n+1) for n in G.nodes()}
    nx.draw(G,pos=pos)
    nx.draw_networkx_labels(G,pos=pos,labels=labels)

    N = G.number_of_nodes()
    edges = np.array(edges,dtype=int)
    A = sprs.csc_matrix((np.ones((edges.shape[0],)),(edges[:,0],edges[:,1])),dtype=float,shape=(N,N))
    A += A.T

    print("recreate results from figure 3 in [1]")
    print("       |\t(i,j)=(1,2)\t(i,j)=(1,3)\t(i,j)=(2,3)")
    print("---------------------------------------------------------")
    for m in range(3):
        T = gtom(A,m)
        print(" m = %d |\t%f\t%f\t%f" %(m,T[0,1],T[0,2],T[1,2]))

    pl.show()

