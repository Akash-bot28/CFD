#!/bin/bash

g++ mesh.cpp -o mesh.o

./mesh.o

python3 plot_beta.py
python3 plot_mesh_spacing.py
python3 plot_mesh.py

