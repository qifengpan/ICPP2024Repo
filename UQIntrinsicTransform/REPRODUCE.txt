This is an description about how to reproduce the data related to the submitted paper of ICPP 2024 for

Improving efficiency of the Monte Carlo method via a code intrinsic framework

This artifact is a detailed implementation of algorithms and integration scheme introduced in section 2 

Requirement:
ROSE compiler https://github.com/rose-compiler/rose
Expected input:
1) The target DC source file, the target language is modern Fortran.
2) Text file that specify the uncertainty source. It should follow the format: Function Name : Variable Name. Indicating the Variable in the function is uncertain.
Expected Outcome:
1) UPP information of the given uncertainty source in text form will be printed. 
2) If user enables the automatic transformation, a source code in UIC form will be generated directly via the toolkit


Example workflow:

1. Build the executable through running: make
2. Go to test code directory :./TestCode/ExplicitHeatSolver_Sequential_DC
3. Generate a uncertainty source description file, as an example, the default text file already exist in the folder as UQ.text
4. running the analysis: ../../UQ_search *.f90 UQ.text

After the search, the UPP information will be directly printed on the screen.
