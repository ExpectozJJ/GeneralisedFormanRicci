"""
A class to compute the Generalised Forman-Ricci curvature for a Vietoris Rips complex from a given NetworkX graph.
"""

# Author:
#     Wee JunJie
#     https://github.com/ExpectozJJ
#
# Remarks:
# Many thanks to stephenhky and saibalmars for their packages MoguTDA and saibalmars. 
# Partial code was adapted from MoguTDA to import simplicial complex from NetworkX.
# stephanhky: https://github.com/stephenhky/MoguTDA
# saibalmars: https://github.com/saibalmars/GraphRicciCurvature

from itertools import combinations
from scipy.sparse import dok_matrix
from operator import add
import numpy as np 
import networkx as nx 

def import_simplices(simplices=[]):
    simplices = map(lambda simplex: tuple(sorted(simplex)), simplices)
    face_set = faces(simplices)
    return list(face_set)
    
def faces(simplices):
    faceset = set()
    for simplex in simplices:
        numnodes = len(simplex)
        for r in range(numnodes, 0, -1):
            for face in combinations(simplex, r):
                faceset.add(tuple(sorted(face)))
    return faceset

def n_faces(face_set, n):
    return filter(lambda face: len(face)==n+1, face_set)

def boundary_operator(face_set, i):
    source_simplices = list(n_faces(face_set, i))
    target_simplices = list(n_faces(face_set, i-1))
    #print(source_simplices, target_simplices)

    if len(target_simplices)==0:
        S = dok_matrix((1, len(source_simplices)), dtype=np.float64)
        S[0, 0:len(source_simplices)] = 1
    else:
        source_simplices_dict = {source_simplices[j]: j for j in range(len(source_simplices))}
        target_simplices_dict = {target_simplices[i]: i for i in range(len(target_simplices))}

        S = dok_matrix((len(target_simplices), len(source_simplices)), dtype=np.float64)
        for source_simplex in source_simplices:
            for a in range(len(source_simplex)):
                target_simplex = source_simplex[:a]+source_simplex[(a+1):]
                i = target_simplices_dict[target_simplex]
                j = source_simplices_dict[source_simplex]
                S[i, j] = -1 if a % 2==1 else 1   # S[i, j] = (-1)**a
    
    return S

class GeneralisedFormanRicci:
    
    def __init__(self, G: nx.Graph):
        """A class to compute Generalised Forman-Ricci curvature for a Vietoris Rips Complex from a given NetworkX graph up to a specified p-dimensional simplex.
        
        Parameters
        ----------
        G : NetworkX graph

        """

        self.G = G.copy()
        self.S = import_simplices(map(tuple, list(nx.find_cliques(self.G))))
        self.v = list(n_faces(self.S, 0))
        self.e = list(n_faces(self.S, 1))
        self.tri = list(n_faces(self.S, 2))
        self.tetra = list(n_faces(self.S, 3))

    def compute_forman(self, p):
        """ Compute Generalised Forman Ricci Curvature for simplices up to dimension p 

        Parameters
        -----------
        p : maximum dimension of simplex to compute curvature.

        Output
        -----------
        forman_dict: dict[simplex, forman curvature value]
            A dictionary of Forman Ricci Curvature values for simplices up to dimension p.
            If p = 0 or p = 1 is entered, it will return the forman ricci curvature for a graph network.
        """

        forman_dict = dict()

        if p < 1:
            p = 1

        bop_data = []
        for i in range(1, p+1):
            bop_data.append(boundary_operator(self.S, i).toarray())

        for simplex in self.S:
            l = len(simplex)
            if l <= p+1 and l >= 2:
                target_simplices = list(n_faces(self.S, l-1))
                target_simplices_dict = {target_simplices[i]: i for i in range(len(target_simplices))}
                m = target_simplices_dict[simplex]
                #print(s, m)
                
                #print(simplex, l)
                b1_ = np.matmul(bop_data[l-1], np.transpose(bop_data[l-1]))-np.diag(np.diag(np.matmul(bop_data[l-1], np.transpose(bop_data[l-1]))))
                b2_ = np.matmul(np.transpose(bop_data[l-2]), bop_data[l-2])-np.diag(np.diag(np.matmul(np.transpose(bop_data[l-2]), bop_data[l-2])))
                p1 = b1_[m]
                p2 = b2_[m]

                forman_dict[l-1][simplex] = sum(bop_data[l-1][m]!=0) + len(simplex) - sum(np.logical_xor(p1,p2))
                #print(simplex, forman_dict[simplex])
        
        for v in self.G.nodes():
            cur_sum = 0
            if self.G.degree(v) > 0:
                for n in self.G.neighbor(v):
                    try:
                        cur_sum += forman_dict[1][(v,n)]
                    except:
                        cur_sum += forman_dict[1][(n,v)]
                forman_dict[0][(v,)] = cur_sum/self.G.degree(v)
            else:
                forman_dict[0][(v,)] = 0
        
        return forman_dict