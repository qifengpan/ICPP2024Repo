This code is originially from https://github.com/fabien-dournac/parallel-heat3d-f90

1. To reproduce the result on Intel 
=============================================
1.1 compile the codes
To modify the mesh size of the code:
1) sequential code: modify size_x ,size_y, size_z in explicitSeq.f90 
2) parallelized code: modify size_x ,size_y, size_z in param

To compile the code: make 
The cmake and openmpi are required, the version is recommend to use openmpi/4.1.2-gcc-11.2.0 to reproduce the results.
successful build will produce both sequential and parallelised executables
=============================================
1.2 measure the performance of the code
1) sequential code
to draw the roofline model of the kernels, please install the intel advisor(2023)
use following command:
advixe-cl --collect=roofline --project-dir=./roofline_results -- ./explicitSeq

the roofline model of all loops and kernels can be read in using GUI of advisor to open the result
2) parallelized code
The elapsed time will be read directly from the output. 

2. To reproduce the result on NEC machine
=============================================
2.1 compile the code
replace the GfORTRAN wit NFORT flag in Make file and then run "make". The current executable will be able to be conducted 
on a NEC machine. On NEC machine, only sequential code will be generated.
2.2 measure the coresponding performance of the kernel
first define the enviroment variable: export VE_PERF_MODE=VECTOR-MEM
Execute the executable: ./explicitSeq
Use ftrace to read the result: ftrace -f ftrace.out

3. Experiment workflow to get the results in the paper:
=================================================
3.1 Sequential code
1) Set the mesh size (see settings in the paper and 1.1)
2) Compiler the code according to instruction in 1.1 on Intel Processor
3) Measure the performance of code on Intel chip via Advisor (see 1.2)
4) Compiler the code on NEC machine (see 2.1), measure performance via ftrace (see 2.2)
5) collect data and change the number of samples on the code do step 2

3.2 Parallelized code
1) Set the mesh size (see settings in the paperi and 1.1)
2) Compiler the code according to instruction in 1.1 on Intel Processor
3) run the code with different number of processors(start from 1 processor)
4) collect the exectution time on different processors and draw the strong scalability curve for given number of samples


