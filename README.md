# General  finite temperature auxiliary field code #

## General information ##
Last version of the code is in Prog_8

Last version of the analysis files is in Analysis_8

1)    Implemented  optimizations made by Johannes.  The code is thereby not necessarily faster !  You  may want to  try the code with the sameoriginal routines upgrade_FFA.f90  Operator_FFA.f90

2)    In the file parameter you can read a variable  CPU_MAX in the namespace &VAR_QMC. If present the code will stop after  the given time. i.e. will  carry out the number of bins you request.  If not present the code will just rect as usual and carry out the requested number of bins. (Implemented by Martin  Bercx). 

3)    In this version,  NWRAP does not have to be a multiple of Ltrot.  This greatly facilitates temperature scans. 

4)    I have included the Condo model on the honeycomb lattice with z-frustrating interactions.

## PREREQUISITES ##

Libraries: Lapack and Blas

Compiler: gfortran  or ifort 


## CONFIGURATION FOR COMPILATION ##


## TESTING ##

## USAGE ##

### INPUT  ###

### OUTPUT ###

### Questions ###

How do we organize the documentation? 

## TODO DOC ##


Johannes Hofmann: documentation for the  SPT model. 

Toshihiro Sato, Martin Hohenadler:  documentation for the Ising model. 

Fakher Assaad:   documentation for the Hubbard model. 

Martin Bercx and Fakher Assaad: general description for the package. 

## TODO CODE ##
cmake

