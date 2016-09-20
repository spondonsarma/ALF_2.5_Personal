export DIR=`pwd`
export f90="gfortran"
export FL="-c -O2 -ftree-vectorize -march=broadwell"
export Libs=${DIR}"/Libraries/"
export LIB_BLAS_LAPACK="-llapack -lblas"
