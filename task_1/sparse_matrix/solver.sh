#!/bin/bash

#icpx sparse_phi.cpp -o sparse_phi
#icpx sparse_psi.cpp -o sparse_psi

./sparse_phi.o
./sparse_psi.o

python3 plot_mesh_phi.py