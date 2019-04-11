#################################
## (c) 2015 Benjamin Maier
#################################
# computes the "Generalized topological overlap" measure as described in [1]
#
# [1] The Generalized Topological Overlap Matrix For
# Detecting Modules in Gene Networks
# by Andy M. Yip, Steve Horvath
# https://labs.genetics.ucla.edu/horvath/GTOM/old/GTOM_tech_report.pdf
from __future__ import print_function


import scipy.sparse as sprs
from numpy import *
import time

def gtom_verbose(A,numSteps,indices=[]):
    return gtom(A,numSteps,indices,verbose=True)

def gtom(A,numSteps,indices=[],verbose=False,callback=None):
    
    #construct matrix B, which encapsulates all neighbors reachable within a path length
    #of numSteps
    number_of_nodes = A.shape[0]
    matrix_shape = A.shape

    #get (NxN) identity
    I = sprs.csr_matrix(sprs.eye(number_of_nodes))

    if numSteps==0:
        return A + I
    else:
        numSteps -= 1

    ###############
    if verbose:
        print("constructing S")
        start = time.time()

    if callback is not None:
        callback("constructing S matrix (this may take a while) ...")
    ###############

    #Get path matrix (which nodes are reachable within numSteps+1 steps)
    S = A
    for m in range(numSteps):
        S = A + S.dot(A)

    ###############
    if verbose:
        end = time.time()
        print("time needed:", end-start)
    ##############


    #get nonzero entries of path matrix
    row,col = S.nonzero()
    no_of_nonzero = len(row)
    del S



    ##############
    if verbose:
        print("construct B")
        print("number of nonzero entries of B:", no_of_nonzero)
        start = time.time()

    if callback is not None:
        callback("constructing B (this may take a while) ...")
    ##############

    #construct B Matrix
    B_data = ones((no_of_nonzero,),dtype=uint32)
    diagonal_entries = nonzero(row==col)[0]
    B_data[diagonal_entries] = 0.
    B = sprs.csr_matrix((B_data,(row,col)),shape=matrix_shape)

    
    ##############
    if verbose:
        end = time.time()
        print("time needed:", end-start)
        print("construct B2")
        start = time.time()
    if callback is not None:
        callback("constructing B^2 ...")
    #############

    #compute the number of reachable nodes by computing B^2
    #added functionality to compute the GTOm-index only for certain nodes
    if len(indices)>0:
        B2_ = B[indices,:].dot(B)
        row,col = B2_.nonzero()
        row = array([ indices[r] for r in row ])
        B2 = sprs.csr_matrix((B2_.data,(row,col)),shape=matrix_shape)
        self_neighborhood = B.sum(axis=1).A1
    else:
        B2 = B.dot(B)
        self_neighborhood = B.sum(axis=1).A1

    #get pairs reachable within 2*(numsteps+1) steps
    row,col = B2.nonzero()

    ##############
    if verbose:
        print("number of nonzero entries of B2:", len(row))
        end = time.time()
        print("time needed:", end-start)
        print("get relevant indices")
        start = time.time()

    if callback is not None:
        callback("getting relevant pairs ...")
    ##############


    #get data for the numerator matrix as described in the paper
    numerator_matrix = B2+A+I
    del I

    #get pairs from the numerator matrix (this is the pairs of nodes for which
    #data is available)
    row,col = numerator_matrix.nonzero()
    numerator_data = numerator_matrix.data
    no_of_nonzero = len(row)


    ##############
    if verbose:
        print("number of FLOs:, ", len(row), ";   expected time:",\
              3.3344888285918445e-06*len(row),"s (",3.3344888285918445e-06*len(row)/60.,"min; or ",\
              3.3344888285918445e-06*len(row)/3600., "h)") 

        end = time.time()
        print("time needed:", end-start)
        print("number of nonzero elements:",no_of_nonzero)
        print("calculating denominator matrix")
        start = time.time()
        
    if callback is not None:
        callback("calculating denominator matrix ...")
    ##############

    #compute the denominator matrix (for element-wise division)
    one = ones((no_of_nonzero,),dtype=float32)    
    denominator_data = one + minimum(self_neighborhood[row],self_neighborhood[col]) - A[row,col].A1

    ##############
    if verbose:
        end = time.time()
        print("time needed:", end-start)
        print("calculating GTOm matrix")
        start = time.time()

    if callback is not None:
        callback("calculating GTOM affinity matrix ...")
    ##############

    #free some memory
    del B2
    del one

    #compute final data
    GTOm_data = numerator_data / denominator_data
    GTOm = sprs.csr_matrix((GTOm_data,(row,col)),shape=matrix_shape)

    ##############
    if verbose:
        end = time.time()
        print("time needed for GTOm matrix:", end-start)
    ##############

    return GTOm




if __name__=="__main__":
    import pylab as pl
    import networkx as nx

    G = nx.Graph()
    G.add_edges_from([(0,1),(1,2),(0,3),(0,4),(0,5),(0,7),
                      (1,3),(1,4),(1,6),(1,8),(1,9),(1,10),
                      (5,6),(7,8)])

    print("nodes:", G.nodes())
    pos = nx.spring_layout(G)
    labels = { n:str(n+1) for n in G.nodes()}
    nx.draw(G,pos=pos)
    nx.draw_networkx_labels(G,pos=pos,labels=labels)
    


    N = G.number_of_nodes()
    A = nx.to_scipy_sparse_matrix(G)

    print("recreate results from figure 3 in [1]")
    print("       |\t(i,j)=(1,2)\t(i,j)=(1,3)\t(i,j)=(2,3)")
    print("---------------------------------------------------------")
    for m in range(3):
        T = gtom(A,m)
        print(" m = %d |\t%f\t%f\t%f" %(m,T[0,1],T[0,2],T[1,2]))

    print() 
    print("computed only for nodes 1 and 2:")
    print("       |\t(i,j)=(1,2)\t(i,j)=(1,3)\t(i,j)=(2,3)")
    print("---------------------------------------------------------")
    for m in range(3):
        T = gtom(A,m,indices=[0,1])
        print(" m = %d |\t%f\t%f\t%f" %(m,T[0,1],T[0,2],T[1,2]))


    #print T[0,1], T[0,2], T[1,2]
    #T = gtom(A,m)
    #print T[10,0], T[10,2], T[10,8]

    pl.show()

    #G2 = nx.Graph()
    #G2.add_edges_from([])
    #pl.show()
