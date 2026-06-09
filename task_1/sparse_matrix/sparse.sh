#!/bin/bash

icpx sparse_3.cpp -o sparse_3

./sparse_3

python3 plot_solver.py 
python3 plot_mesh_phi.py $id


