export DIR=`pwd`
export f90="gfortran"
export F90OPTFLAGS="-O3 -march=broadwell"
export FL="-c ${F90OPTFLAGS}"
export Libs=${DIR}"/Libraries/"
export LIB_BLAS_LAPACK="-llapack -lblas"
#export LIB_BLAS_LAPACK="-Wl,--start-group /opt/intel/compilers_and_libraries_2017.0.098/linux/mkl/lib/intel64_lin/libmkl_gf_lp64.a /opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64_lin/libmkl_core.a /opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64_lin/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl"
