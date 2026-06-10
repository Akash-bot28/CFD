#!/bin/bash

icpx sparse_2.cpp -o sparse_2
icpx sparse_psi.cpp -o sparse_psi

./sparse_2
./sparse_psi

#python3 plot_solver.py 

python3 plot_mesh_phi.py


