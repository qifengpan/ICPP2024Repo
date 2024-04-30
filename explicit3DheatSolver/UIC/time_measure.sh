#!/bin/bash
start=$(date +%s%N) 

mpirun -np 32 explicitPar

end=$(date +%s%N)
time=$((end-start))
echo "time used:$time seconds"

