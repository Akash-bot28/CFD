#!/bin/bash

#icpx mesh.cpp -o mesh.o
#icpx sparse_phi.cpp -o sparse_phi.o
#icpx sparse_psi.cpp -o sparse_psi.o

#./mesh
./sparse_phi.o
./sparse_psi.o

python3 plot_solver_convergence.py