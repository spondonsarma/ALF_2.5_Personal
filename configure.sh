export DIR=`pwd`
#export f90=$mpif90
export f90="gfortran"
export FL="-c  -w  -O3"
export Libs=${DIR}"/Libraries/"
export LIB_BLAS_LAPACK="-llapack -lblas"
