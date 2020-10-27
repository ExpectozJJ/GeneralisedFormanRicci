[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

# GeneralisedFormanRicci
This code computes the Forman Ricci Curvature for simplicial complex generated from a given point cloud data. The implementation is based on the combinatorial definition of Forman Ricci curvature defined by Robin Forman. This implementation generalises beyond the simplified version implemented in saibalmars/GraphRicciCurvature github.

Many thanks to stephenhky and saibalmars for their packages MoguTDA and GraphRicciCurvature respectively. 
Partial code was modified from MoguTDA for the computation of the boundary matrices. 

## Installation

`pip install GeneralisedFormanRicci`

Upgrading via `pip install --upgrade GeneralisedFormanRicci`

## Package Requirement

* [NetworkX](https://github.com/networkx/networkx) >= 2.0 (Based Graph library)
* [GUDHI](https://github.com/GUDHI) (Simplicial Complex Library)
* [NumPy](https://github.com/numpy/numpy)

## Simple Example

```
from GeneralisedFormanRicci.frc import GeneralisedFormanRicci

data = [[0.8, 2.6], [0.2, 1.0], [0.9, 0.5], [2.7, 1.8], [1.7, 0.5], [2.5, 2.5], [2.4, 1.0], [0.6, 0.9], [0.4, 2.2]]
for f in [0, 0.5, 1, 2, 3]:
    sc = GeneralisedFormanRicci(data, method = "rips", epsilon = f)
    sc.compute_forman()
```

## References
* MoguTDA: https://github.com/stephenhky/MoguTDA
* GraphRicciCurvature: https://github.com/saibalmars/GraphRicciCurvature
* Forman, R. (2003). Bochner's method for cell complexes and combinatorial Ricci curvature. Discrete and Computational Geometry, 29(3), 323-374.
* Forman, R. (1999). Combinatorial Differential Topology and Geometry. New Perspectives in Algebraic Combinatorics, 38, 177.

## Cite 
If you use this code in your research, please considering cite our paper:
