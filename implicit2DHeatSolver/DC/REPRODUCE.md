This code is originiall from: https://github.com/omersan/5.03.Heat2D/tree/master

1. To reprocude the Result on Intel processor:
================================================
1.1 compile the code
using compiling flags:
gfortran -O3 -g -march='cascadelake' -mprefer-vector-width=512 -funroll-loops *.f90 -o heatTrans
================================================
1.2 measure the coresponding performance of the kernel
to draw the roofline model of the kernels, please install the intel advisor(2023)
use following command:
advixe-cl --collect=roofline --project-dir=./roofline_results -- ./heatTrans.f90

the roofline model of all loops and kernels can be read in using GUI of advisor to open the result

2. To reproduce the results on NEC vector machine:
=================================================
2.1 compile the code
No code modification is needed, directly use the code same as used above.
using compiling flags:
nfort -Wall -O3 -ftrace -floop-unroll-complete-nest=3 *.f90 -o heatTrans_Au
================================================
2.2 measure the coresponding performance of the kernel
first define the enviroment variable: export VE_PERF_MODE=VECTOR-MEM
Execute the executable: ./heatTrans_Au
Use ftrace to read the result: ftrace -f ftrace.out

3. Experiment workflow to get the results in the paper:
=================================================
1) Compiler the code according to instruction in 1.1 on Intel Processor
2) Measure the performance of code on Intel chip via Advisor (see 1.2)
3) Compiler the code on NEC machine (see 2.1), measure performance via ftrace (see 2.2)
4) collect data 
