
from gtom import gtom
from gtom.gtom_calculator import GTOMCalc
import scipy.sparse as sprs
import numpy as np

import pylab as pl
import networkx as nx

edges = [
            (0,1),(1,3),(1,4),(1,5),(1,11),
            (3,6),(4,8),(5,9),(10,11),
            (7,2),(2,6),(2,8),(2,9),(2,10),
        ]

G = nx.Graph()
G.add_edges_from(edges)
print("nodes:", G.nodes())
pos = nx.spring_layout(G)
labels = { n:str(n) for n in G.nodes()}
nx.draw(G,pos=pos)
nx.draw_networkx_labels(G,pos=pos,labels=labels)



N = G.number_of_nodes()

edges = np.array(edges,dtype=int)
A = sprs.csc_matrix((np.ones((edges.shape[0],)),(edges[:,0],edges[:,1])),dtype=float,shape=(N,N))
A += A.T
print(A.todense())

print("recreate results from figure 3 in [1]")
print("       |\t(i,j)=(1,2)\t(i,j)=(1,3)\t(i,j)=(2,3)")
print("---------------------------------------------------------")
for m in range(3):
    T = gtom(A,m)
    print(" m = %d |\t%f\t%f\t%f" %(m,T[1,2],T[1,3],T[2,3]))

print() 
print("computed only for nodes 1 and 2:")
print("       |\t(i,j)=(1,2)\t(i,j)=(1,3)\t(i,j)=(2,3)")
print("---------------------------------------------------------")
for m in range(3):
    gtom_calc = GTOMCalc(A,m)
    T = gtom_calc.compute_for_indices(indices=[1,2])
    print(" m = %d |\t%f\t%f\t%f" %(m,T[1,2],T[1,3],T[2,3]))


pl.show()
