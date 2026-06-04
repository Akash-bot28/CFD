#!/bin/bash

#icpx mesh.cpp -o mesh
./mesh

python3 plot_d.py
python3 plot_mesh.py