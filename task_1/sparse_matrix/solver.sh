#!/bin/bash

#icpx sparse_phi.cpp -o sparse_phi
#icpx sparse_psi.cpp -o sparse_psi

./sparse_phi
./sparse_psi

python3 plot_mesh_phi.py