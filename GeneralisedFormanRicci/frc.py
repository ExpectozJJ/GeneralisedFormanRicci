"""
A class to compute the Generalised Forman-Ricci curvature for a Simplicial Complex from a given point cloud data.
"""

# Author:
#     Wee JunJie
#     https://github.com/ExpectozJJ
#
# Remarks:
# Many thanks to stephenhky and saibalmars for their packages MoguTDA and saibalmars. 
# Partial code was adapted from MoguTDA to import simplicial complex.
# stephanhky: https://github.com/stephenhky/MoguTDA
# saibalmars: https://github.com/saibalmars/GraphRicciCurvature

from itertools import combinations
from scipy.sparse import dok_matrix
from scipy.spatial import Delaunay
from operator import add, itemgetter
import numpy as np 
import networkx as nx 
import gudhi

from collections import defaultdict


def gen_graph(edges, points, labels):
    G = nx.Graph()
    for i in range(len(points)):
        G.add_node(i, coords = points[i])
        if len(labels.keys()) > 0:
            for l in labels.keys():
                G.nodes[i][l] = labels[l][i]

    for e in edges:
        G.add_edge(e[0], e[1])

    return G
    
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
                S[i, j] = -1 if a % 2==1 else 1
    
    return S


class GeneralisedFormanRicci:
    
    def __init__(self, points, labels=None, epsilon=2.0, method="rips", p = 2):
        """A class to compute Generalised Forman-Ricci curvature for a specified p-dimensional simplex from a simplicial complex generated from point cloud data
        
        Parameters
        ----------
        points: n-dimensional point cloud data
        labels: dictionary of attributes for points
        method: 
            Type of Simplicial Complex to be generated
            -------------------------------------------
            rips: Vietoris Rips Complex
            alpha: Alpha Complex
        epsilon: Diameter

        """

        self.pts = np.array(points)
        self.labels = {'coords': self.pts} if labels == None else labels
        self.method = method
        self.epsilon = epsilon
        self.p = p

        if self.method == "rips":
            try:
                self.S = self.construct_rips(self.p)
                print("Rips Complex Constructed.")
            except:
                raise('epsilon not defined for rips method.')
        elif self.method == "alpha":
            self.S = self.construct_alpha(self.p)
        else:
            raise("Unknown method specified.")


        """
        Compute and Store Hodge Laplacian up to dimension p. 
        """
        self.laplacian = []
        self.laplacian.append(np.matmul(boundary_operator(self.S, 1).toarray(), np.transpose(boundary_operator(self.S, 1).toarray())))
        for i in range(1, p+2):
            b1 = boundary_operator(self.S, i+1).toarray()
            b2 = boundary_operator(self.S, i).toarray()
            b1_ = np.matmul(b1,np.transpose(b1))
            b2_ = np.matmul(np.transpose(b2), b2)
            self.laplacian.append(b1_+b2_)

        # Constructing Networkx Graph might take awhile.
        #self.G = gen_graph(list(n_faces(self.S, 1)), self.pts, self.labels)

    def construct_alpha(self, p):
        alpha_complex = gudhi.AlphaComplex(self.pts)
        alpha = alpha_complex.create_simplex_tree()
        val = alpha.get_filtration()
        simplices = set()
        for v in val:
            if len(v[0]) <= p+1 and np.sqrt(v[1])*2 <= self.epsilon: #circumradius must be converted to diameter and within the filtration parameter
                simplices.add(tuple(v[0]))

        return simplices

    def construct_rips(self, p):
        rips_complex = gudhi.RipsComplex(points = self.pts, max_edge_length = self.epsilon)
        simplex_tree = rips_complex.create_simplex_tree(max_dimension = p+1)
        val = simplex_tree.get_filtration()
        simplices = set()
        for v in val:
            if len(v[0]) <= p+1:
                simplices.add(tuple(v[0]))

        return simplices

    def _compute_forman(self, simplex):
        """
        Lookup the values in Hodge Laplacian and Compute the Forman Ricci Curvature
        """
        l = len(simplex)
        if l <= self.p+1 and l >= 2:
            target_simplices = list(n_faces(self.S, l-1))
            target_simplices_dict = {target_simplices[i]: i for i in range(len(target_simplices))}
            m = target_simplices_dict[simplex]
            
            """
            Minimise usage of np.matmul and np.diag as time takes longer for large simplicial complex.
            """
            mat = np.delete(self.laplacian[l-1][m], m) 
            val = self.laplacian[l-1][m][m] - sum(np.abs(mat))
            
            return val

    def compute_forman(self):
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

        forman_dict = defaultdict(dict)

        p = self.p
        if p < 1:
            p = 1

        
        simplices = list(self.S)
        
        for i in range(len(simplices)):
            l = len(simplices[i])
            if l <= p+1 and l >= 2:
                forman_dict[l-1][simplices[i]] = self._compute_forman(simplices[i])
        
        target_simplices = list(n_faces(self.S, 0))
        cnt = np.zeros(len(target_simplices))
        vals = np.zeros(len(target_simplices))
        nbrs = list(n_faces(self.S, 1))
        for vertex in nbrs:
            for v in vertex:
                vals[v] += forman_dict[1][vertex]
                cnt[v] += 1
        for i in range(len(vals)):
            if cnt[i] > 0:
                forman_dict[0][(i,)] = vals[i]/cnt[i]
            else:
                forman_dict[0][(i,)] = 0
        
        return forman_dict