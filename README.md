Ising_Hypercube
===============

A Monte Carlo implementation of the Ising Hamiltonian on a hypercube of arbitrary dimension

| D | Tc |
|----|----|
|2 | 2.27|
|3 | 4.51|
|6 | 10.83|
|7 | 12.87|

An explanation of the cell-counting:
https://en.wikipedia.org/wiki/Hypercube


Percolation
-----------
http://arxiv.org/abs/cond-mat/0101295

Number of cells in general D
-----------

| D | N_1 | N_2 |  N_3 |
|----|----|---|---|
|1 | 1 | 0 | 0 |
|2 | 2 | 1 | 0 |
|3 | 3 | 3 | 1 |
|4 | 4 | 6 | 4 |
|D | Dc1 | Dc2 | Dc3 |
