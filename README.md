
## General information ##
This version of the **A**lgorithms for **L**attice **F**ermions package provides a general code for the finite temperature auxiliary field Quantum Monte Carlo algorithm.       The code  is engineered to  be able simulate any model that can be written in terms of  sums of single body operators, of squares of single body operators and single body operators coupled to an Ising field with  given dynamics. We  provide predefined types that allow  the user to specify the model, the  Bravais lattice  as well as equal time and time displaced observables.     The code supports an MPI implementation.   Examples such as the Hubbard model on the Honeycomb lattice  as well as the Hubbard model  on the square lattice coupled to a transverse Ising field are  provided and discussed in the [documentation](https://git.physik.uni-wuerzburg.de/fassaad/General_QMCT_code/raw/documentation_new/Documentation/doc.pdf).  

The Hamiltonians we can consider read:
![The Hamiltonian0](https://git.physik.uni-wuerzburg.de/fassaad/General_QMCT_code/raw/master/Images/Hamiltonian0.png)
where
![The Hamiltonian1](https://git.physik.uni-wuerzburg.de/fassaad/General_QMCT_code/raw/master/Images/Hamiltonian1.png)
where Z denotes an Ising spin variable with predefined dynamics. If your model can be written in this, form then it will be amenable to the ALF. 
## PREREQUISITES ##

Libraries: Lapack and Blas

Compiler: gfortran  or ifort 


## CONFIGURATION FOR COMPILATION ##
Edit **set_env.sh**   to provide compiler information and flags as well as  location of Lapack and Blas.

## FILES AND DIRECTORIES ##

**Libraries**    Libraries. Once that the environment is set in the file set_env.sh  the Libraries can be compiled with the **make** command. 

**Prog**   Main program and subroutines.   The command **make Examples** will  generate the **Examples.out** executable. 

**Analysis** Analysis programs. The command **make** will generate three executables required  for the error analysis  of scalar, equal time and time displaced observables. 

**Start**   This directory contain the files required to start a run. In particular it contains the parameter file   that specifies the model the lattice and various   parameters for the Monte Carlo run and  error analysis. 

**Examples** This directory provides a set of short example runs.  

**Documentation**  We have included in the file ALF-05.pdf an extensive documentation, 
 

## TESTING ##

We have about 30 tests that test various parts of the program in the folder testsuite.
As testing framework we employ CTest.
From the subfolder testsuite the tests can be run as follows
- mkdir tests
- cd tests
- cmake ..
- make
- make test


## LICENSE ##
The various works that make up the ALF project are placed under licenses that put
a strong emphasis on the attribution of the original authors and the sharing of the contained knowledge.
To that end we have placed the ALF source code under the GPL version 3 license (see license.GPL and license.additional)
and took the liberty as per GPLv3 section 7 to include additional terms that deal with the attribution
of the original authors(see license.additional).
The Documentation of the ALF project by the ALF contributors is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License (see Documentation/license.CCBYSA)
We mention that we link against parts of lapack which licensed under a BSD license(see license.BSD).

