# Fortran_DNAD
Implement DNAD (Automatic Differentiation) for Fortran Code

This code uses combines 3 items already developed in the literature:
1. John Newman's `BAND` subroutine (found in Appendix C page 619 of his book Electrochemical Systems 3rd Ed. (2004))
2. `fillmat` found in Chapter 7 (page 129) of J. Deliang Yang's these (Columbia University 1997 [thesis](https://clio.columbia.edu/catalog/1987854?counter=1))
3. And the dual number automatic differentiation Fortran module `dnadmod` (originally published as "DNAD, a Simple Tool for Automatic
Differentiation of Fortran Codes Using Dual Numbers" by Wenbin Yu and Maxwell Blair; and can be found [here](https://www.sciencedirect.com/science/article/pii/S0010465513000027?casa_token=MpXIh34txb0AAAAA:vf9mYSrbAU3VNKE9MYdLnQkd2OpTSa2AW0D5sN9FNbCI9fkhPZw-UXEcbR_4-CYoKAwEXgXmivA))
