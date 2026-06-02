#!/bin/bash


#icpx mesh.cpp -o mesh.o

./mesh.o

python3 plot_beta.py $id
python3 plot_d.py $id
python3 plot_mesh.py $id

