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

def _compute_forman(faceset, simplex):
    """ Computes the forman curvature for simplex input with a given simplicial complex. 

    Inputs
    -----------
    faceset: simplicial complex 
    simplex: i-dimensional simplex 

    Output: 
    -----------
    Forman Ricci Curvature for i-dimensional simplex
    """

    i = len(simplex)
    if i >= 2:
        #source_simplices = list(n_faces(self.S, i))
        target_simplices = list(n_faces(faceset, i-1))

        #source_simplices_dict = {source_simplices[j]: j for j in range(len(source_simplices))}
        target_simplices_dict = {target_simplices[i]: i for i in range(len(target_simplices))}

        b1 = boundary_operator(faceset, i).toarray()
        b2 = boundary_operator(faceset, i-1).toarray()
        b1_ = np.matmul(b1, np.transpose(b1))-np.diag(np.diag(np.matmul(b1, np.transpose(b1))))
        b2_ = np.matmul(np.transpose(b2), b2)-np.diag(np.diag(np.matmul(np.transpose(b2), b2)))
        p1 = b1_[target_simplices_dict[simplex]]
        p2 = b2_[target_simplices_dict[simplex]]

        return sum(b1[target_simplices_dict[simplex]]!=0) + len(simplex) - sum(np.logical_xor(p1,p2))

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

        if p > 1:
            for i in range(len(self.S)):
                #print(self.S[i], _compute_prl_nbr(self.S, self.S[i]))
                if len(self.S[i])-1 <= p and len(self.S[i])-1 > 0:
                    forman_dict[self.S[i]] = _compute_forman(self.S, self.S[i])
        else:
            for (v1, v2) in self.G.edges():
                forman_dict[(v1,v2)] = _compute_forman(self.S, (v1,v2))
        
        for v in self.G.nodes():
            cur_sum = 0
            if self.G.degree(v) > 0:
                for n in self.G.neighbor(v):
                    try:
                        cur_sum += forman_dict[(v,n)]
                    except:
                        cur_sum += forman_dict[(n,v)]
                forman_dict[(v,)] = cur_sum/self.G.degree(v)
            else:
                forman_dict[(v,)] = 0
        
        return forman_dict