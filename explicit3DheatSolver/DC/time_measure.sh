#!/bin/bash
start=$(date +%s.%3N) 

for i in {1..128}
do
#    mpirun -np 32 explicitPar 
     ./explicitSeq
done

end_t=$(date +%s.%3N)
#time=$((end-start))
time = $(echo "scale=3; $end_t - $start" | bc)
echo "time used:$time seconds"
#timediff $start $end
