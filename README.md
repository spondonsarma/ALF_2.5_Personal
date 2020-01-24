# ALF #
[![pipeline status](https://git.physik.uni-wuerzburg.de/fassaad/General_QMCT_code/badges/master/pipeline.svg)](https://git.physik.uni-wuerzburg.de/fassaad/General_QMCT_code/commits/master)
[![coverage report](https://git.physik.uni-wuerzburg.de/fassaad/General_QMCT_code/badges/master/coverage.svg)](https://git.physik.uni-wuerzburg.de/fassaad/General_QMCT_code/commits/master)
## General information ##
This version of the **A**lgorithms for **L**attice **F**ermions package provides a general code for the finite temperature auxiliary field Quantum Monte Carlo algorithm.       The code  is engineered to  be able simulate any model that can be written in terms of  sums of single body operators, of squares of single body operators and single body operators coupled to an Ising field with  given dynamics. We  provide predefined types that allow  the user to specify the model, the  Bravais lattice  as well as equal time and time displaced observables.     The code supports an MPI implementation.   Examples such as the Hubbard model on the Honeycomb lattice  as well as the Hubbard model  on the square lattice coupled to a transverse Ising field are  provided and discussed in the [documentation](https://git.physik.uni-wuerzburg.de/fassaad/General_QMCT_code/blob/master/Documentation/doc.pdf) and can be tested using ALF's [Test Suite](https://git.physik.uni-wuerzburg.de/fassaad/Testsuite_General_QMCT_code).

The Hamiltonians we can consider read:
```math
\hat{\mathcal{H}}=\hat{\mathcal{H}}_{T}+\hat{\mathcal{H}}_{V} +  \hat{\mathcal{H}}_{I} +   \hat{\mathcal{H}}_{0,I}
```
where
```math
\begin{aligned}
\hat{\mathcal{H}}_{T}
&=
\sum\limits_{k=1}^{M_T}
\sum\limits_{\sigma=1}^{N_{\mathrm{col}}}
\sum\limits_{s=1}^{N_{\mathrm{fl}}}
\sum\limits_{x,y}^{N_{\mathrm{dim}}}
\hat{c}^{\dagger}_{x \sigma   s}T_{xy}^{(k s)} \hat{c}^{\phantom\dagger}_{y \sigma s}\\
\hat{\mathcal{H}}_{V}
&=
\sum\limits_{k=1}^{M_V}U_{k}
\left\{
\sum\limits_{\sigma=1}^{N_{\mathrm{col}}}
\sum\limits_{s=1}^{N_{\mathrm{fl}}}
\left[
\left(
\sum\limits_{x,y}^{N_{\mathrm{dim}}}
\hat{c}^{\dagger}_{x \sigma s}V_{xy}^{(k s)}\hat{c}^{\phantom\dagger}_{y \sigma s}
\right)
+\alpha_{k s}
\right]
\right\}^{2} \\
\hat{\mathcal{H}}_{I}
& = 
\sum\limits_{k=1}^{M_I} \hat{Z}_{k}
\left(
\sum\limits_{\sigma=1}^{N_{\mathrm{col}}}
\sum\limits_{s=1}^{N_{\mathrm{fl}}}
\sum\limits_{x,y}^{N_{\mathrm{dim}}}
\hat{c}^{\dagger}_{x \sigma s} I_{xy}^{(k s)}\hat{c}^{\phantom\dagger}_{y \sigma s}
\right) 
\;.
\end{aligned}

```

Here Z denotes an Ising spin variable with predefined dynamics. If your model can be written in this form then it will be amenable to the ALF. 

## Doxygen ##

You can find here [Doxygen](https://pawn.physik.uni-wuerzburg.de/~assaad/Doxygen_Docu/ALF/html/index.html)  formatted documentation. (Work in progress)

## PREREQUISITES ##

Libraries: Lapack and Blas

Compiler: gfortran or ifort 


## CONFIGURATION FOR COMPILATION ##

**configureHPC.sh**  It is recommended to use this script to set the environment variables. Type ./configureHPC.sh to  browse through a list of options.

Once you have run the configuration script, run **make** to compile the whole package. Alternatively, run the Makefiles separately by first changing to directory Libraries and then Analysis. In the Prog directory then type make examples.   The other programs are being updated to comply with the new version of the code.

## FILES AND DIRECTORIES ##

**Libraries**   Libraries.

**Prog**   Main program and subroutines.

**Analysis**   Analysis programs. 

**Scripts_and_Parameters_files**  Suite of helper scripts and the Start directory which contains the files required to start a run, in particular  the **parameters** file that specifies the model, the lattice, and various parameters for the Monte Carlo run and  error analysis. 

**Documentation**   We have included in the file  [doc.pdf](https://git.physik.uni-wuerzburg.de/fassaad/General_QMCT_code/blob/master/Documentation/doc.pdf)   an extensive documentation. The development of the documentation will take place **only** in the documentation_new branch. In this branch you all have permission to push things so that we can keep on improving  things. 
 

## TESTING ##

We have about 30 tests that test various parts of the program in the folder testsuite.
As testing framework we employ CTest.
From the ALF folder the tests can be run as follows
```bash
source configureHPC.sh Devel serial
gfortran -v
make lib
make ana
make Examples
cd testsuite
cmake -E make_directory tests
cd tests
cmake -G "Unix Makefiles" -DCMAKE_Fortran_FLAGS_RELEASE=${F90OPTFLAGS} \
      -DCMAKE_BUILD_TYPE=RELEASE ..
cmake --build . --target all --config Release
ctest -VV -O log.txt
```


## LICENSE ##
The various works that make up the ALF project are placed under licenses that put
a strong emphasis on the attribution of the original authors and the sharing of the contained knowledge.
To that end we have placed the ALF source code under the GPL version 3 license (see license.GPL and license.additional)
and took the liberty as per GPLv3 section 7 to include additional terms that deal with the attribution
of the original authors(see license.additional).
The Documentation of the ALF project by the ALF contributors is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License (see Documentation/license.CCBYSA)
We mention that we link against parts of lapack which licensed under a BSD license(see license.BSD).

## KNOWN ISSUES ##

Intel suite 2017: We have detected a bug in both 2017 versions of Intel's MPI implementation (mpi.intel/2017(default) 
and mpi.intel/2017.2) if used in combination with the parallel (threaded) MKL library. The advice is to 
use either the 2016 suite (intel/16.0 (compiler), mpi.intel/5.1 and mkl/11.3 ) or the new 2018 suite 
(intel/18.0 (compiler), mpi.intel/2018 and mkl/2018). We did not detect this issue in either of the environments. 
You should also be aware that, by default, dynamic linking is used. Hence if you use the 2016 or 2018 modules 
at compilations, the bug can reenter if you still load the 2017 versions at runtime. So please adapt your
configureHPC.sh as well as your Jobfiles for the loadleveler accordingly.
Additional note: In the serial version, the bug also seems to be absent. 
If you want to use the 2017 suite, you have to use the serial version of MKL (mkl/2017_s), which means you 
cannot profit from openMP multi-threading. This library is linked statically, hence taking care of this at 
compile time is sufficient and there is no need to adapt the Jobfiles.
WARNING: Even if you do not use parallel tempering actively, we still strongly suggest to take care of 
the above bug as it is extremely hard to estimate hidden influences and correlations of this memory 
allocation bug in the rest of the program. It is possible the other parts of the algorithm might be 
affected besides the tempering exchange step even so we have absolutely no hint of additionally 
affected sections in ALF.


    

