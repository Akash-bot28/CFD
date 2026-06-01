#!/bin/bash

id=$1

#icpx mesh_${id}.cpp -o mesh_${id}.o

./mesh_${id}.o

python3 plot_beta.py $id
python3 plot_d.py $id
python3 plot_mesh.py $id

