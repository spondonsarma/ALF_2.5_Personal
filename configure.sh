# setting QRREF has the highest priority. Setting nothing selects System lapack for the QR decomposition
PROGRAMMCONFIGURATION="-DQRREF"
DIR=`pwd`
export DIR
f90="gfortran"
export f90
FL="-c -O3 ${PROGRAMMCONFIGURATION}"
export FL
Libs=${DIR}"/Libraries/"
export Libs
LIB_BLAS_LAPACK="-llapack -lblas"
export LIB_BLAS_LAPACK
