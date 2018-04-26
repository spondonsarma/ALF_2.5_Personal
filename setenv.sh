ALF_FC=$mpif90
ALF_FC=gfortran
F90OPTFLAGS="-O3 -Wconversion -fcheck=all"
F90OPTFLAGS="-O3 -Wconversion"
F90USEFULFLAGS="-cpp -std=f2003"
#F90USEFULFLAGS=" -cpp"
LIB_BLAS_LAPACK="-llapack -lblas"

export F90USEFULFLAGS
FL="-c ${F90OPTFLAGS} ${PROGRAMCONFIGURATION}"
export FL
Libs="$DIR/Libraries/"
export Libs
export LIB_BLAS_LAPACK

Libs="$(pwd)/Libraries"
ALF_INC="-I${Libs}/Modules ${INC_HDF5}"
export ALF_LIB="${Libs}/Modules/modules_90.a ${LIB_BLAS_LAPACK} ${Libs}/libqrref/libqrref.a"

export ALF_DIR="$(pwd)"
export ALF_FC

export ALF_FLAGS_PROG="-c ${F90USEFULFLAGS} ${F90OPTFLAGS} ${PROGRAMCONFIGURATION} ${ALF_INC}"
export ALF_FLAGS_QRREF="-c ${F90OPTFLAGS}"
export ALF_FLAGS_ANA="-c -cpp ${F90OPTFLAGS} ${ALF_INC}"
export ALF_FLAGS_MODULES="-c -cpp ${F90OPTFLAGS}"
