export DIR=`pwd`
export f90=$mpif90
export f90="gfortran"
export FL="-c  -w   -r8  -O3   -g  -pg"
export FL="-c  -w   -r8  -O3  -fbounds-check -ftrace=full"
export FL="-c  -w   -r8  -O3 "
export FL="-c  -w   -fdefault-real-8 -O3 "
export FL="-c  -w   -fdefault-real-8  -fbounds-check -g -fbacktrace -ffpe-trap=zero,invalid"
export FL="-c  -w  -O3"
export Libs=${DIR}"/Libraries/"
export LIB_BLAS_LAPACK=${HOME}"/lib_90/LaPack/lapack.a  "${HOME}"/lib_90/Blas/libblas.a"
