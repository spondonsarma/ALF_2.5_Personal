# setting QRREF has the highest priority. Setting nothing selects System lapack for the QR decomposition.
# Setting OLDNAG selects syntax of NAG Versions before NAG Mark 17 (roughly....).
# In addition the NAG library has to be specified in LIB_BLAS_LAPACK.
# -DMPI selects MPI.
PROGRAMMCONFIGURATION="-DQRREF -DSTAB1"
PROGRAMMCONFIGURATION=""
f90="ifort"
f90="/opt/pgi/linux86-64/2016/bin/pgfortran"
export f90
F90OPTFLAGS="-O3 "
F90OPTFLAGS="-O3 -Wconversion  -fcheck=all"
F90OPTFLAGS="-O3 -Minform=inform"
F90USEFULFLAGS="-Mpreprocess -Minform=inform"
export F90USEFULFLAGS
export F90OPTFLAGS
FL="-c ${F90OPTFLAGS} ${PROGRAMMCONFIGURATION}"
export FL
DIR=`pwd`
export DIR
Libs=${DIR}"/Libraries/"
export Libs
LIB_BLAS_LAPACK="-L/opt/pgi/linux86-64/2016/lib -llapack -lblas"
export LIB_BLAS_LAPACK
