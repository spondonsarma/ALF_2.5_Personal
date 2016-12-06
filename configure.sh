export DIR=`pwd`
export f90="gfortran"
export FL="-c -O3"
export Libs=${DIR}"/Libraries/"
export LIB_BLAS_LAPACK="-llapack -lblas"
export Enable_MPI=false

if $Enable_MPI ; then
   echo '#define MPI' > Prog_8/machine
else
   echo '#define noMPI' > Prog_8/machine
fi
