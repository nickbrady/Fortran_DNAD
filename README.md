# Fortran_DNAD
Implement DNAD (Automatic Differentiation) into Newman's `BAND` subroutine

This code combines 3 items already developed in the literature:
1. John Newman's `BAND` subroutine (found in Appendix C page 619 of his book Electrochemical Systems 3rd Ed. (2004))
2. `fillmat` found in Chapter 7 (page 129) of the thesis: Experimental and numerical investigation of mass transfer in electrochemical systems by J. Deliang Yang (Columbia University, 1997 [link](https://clio.columbia.edu/catalog/1987854?counter=1))
3. And the dual number automatic differentiation Fortran module `dnadmod` (originally published as "DNAD, a Simple Tool for Automatic
Differentiation of Fortran Codes Using Dual Numbers" by Wenbin Yu and Maxwell Blair (found [here](https://www.sciencedirect.com/science/article/pii/S0010465513000027?casa_token=MpXIh34txb0AAAAA:vf9mYSrbAU3VNKE9MYdLnQkd2OpTSa2AW0D5sN9FNbCI9fkhPZw-UXEcbR_4-CYoKAwEXgXmivA)); `dnadmod` was originally forked from joddlehod's [repository](https://github.com/joddlehod/dnad))

The use of automatic differentiation significantly simplifies the process of linearizing systems of differential equations without a large performance loss and maintains numerical accuracy to within machine precision.

# Example Problems
Some example problems are included to help the user understand how models can be developed within this framework

## Code Structure
Because `dnadmod` is a long module and does not need to be edited frequently, it is logical to put it at the end of the Fortran files. However, Fortran necessitates that modules appear before subsequent modules that depend on them. 
To be able to keep `dnadmod` at the end of the file and also be able to run the program, the python script "SortFortranModules.py" was written to re-order the modules and then run the programs. 
If one is not interested in using a python script to run the Fortran programs, one can elect to run the files that begin with "ReSort" from the terminal by first compiling the program: `gfortran -fdefault-real-8 _file_name_` and then executing the program: `./a.out`
