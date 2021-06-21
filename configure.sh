#!/bin/sh
STABCONFIGURATION=""
# STABCONFIGURATION="${STABCONFIGURATION} -DQRREF"

export ALF_DIR="$PWD"

set_hdf5_flags()
{
  CC="$1" FC="$2" CXX="$3"
  
  $FC -o get_compiler_version.out get_compiler_version.F90
  compiler_vers=$(./get_compiler_version.out | sed 's/ /_/g')
  
  HDF5_DIR="$ALF_DIR/HDF5/$compiler_vers"
  if [ ! -d "$HDF5_DIR" ]; then
    printf "\e[31mDownloading and installing HDF5 in %s.\e[0m\n" "$HDF5_DIR"
    CC="$CC" FC="$FC" CXX="$CXX" HDF5_DIR="$HDF5_DIR" "$ALF_DIR/HDF5/install_hdf5.sh" || return 1
  fi
  INC_HDF5="-I$HDF5_DIR/include"
  LIB_HDF5="-L$HDF5_DIR/lib $HDF5_DIR/lib/libhdf5hl_fortran.a $HDF5_DIR/lib/libhdf5_hl.a"
  LIB_HDF5="$LIB_HDF5 $HDF5_DIR/lib/libhdf5_fortran.a $HDF5_DIR/lib/libhdf5.a -lz -ldl -lm -Wl,-rpath -Wl,$HDF5_DIR/lib"
}

# default optimization flags for Intel compiler
INTELOPTFLAGS="-cpp -O3 -fp-model fast=2 -xHost -unroll -finline-functions -ipo -ip -heap-arrays 1024 -no-wrap-margin"
INTELOPTFLAGS="-cpp -O3 "
#INTELOPTFLAGS="$INTELOPTFLAGS -traceback"
# uncomment the next line if you want to use additional openmp parallelization
INTELOPTFLAGS="${INTELOPTFLAGS} -parallel -qopenmp"
INTELDEVFLAGS="-warn all -check all -g -traceback"
INTELUSEFULFLAGS="-std08"

# default optimization flags for GNU compiler
GNUOPTFLAGS="-cpp -O3 -ffree-line-length-none -ffast-math"
#GNUOPTFLAGS="-cpp -O0 -ffree-line-length-none"
# uncomment the next line if you want to use additional openmp parallelization
GNUOPTFLAGS="${GNUOPTFLAGS} -fopenmp"
# GNUDEVFLAGS="-Wconversion -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow,denormal"
GNUDEVFLAGS="-Wconversion -Werror -fcheck=all -g -fbacktrace -fmax-errors=10"
GNUUSEFULFLAGS="-std=f2008"

# default optimization flags for PGI compiler
PGIOPTFLAGS="-Mpreprocess -O1"
# uncomment the next line if you want to use additional openmp parallelization
PGIOPTFLAGS="${PGIOPTFLAGS} -mp"
PGIDEVFLAGS="-Minform=inform -g -traceback"
PGIUSEFULFLAGS=""

MACHINE=""
Machinev=0
MODE=""
modev=0
STAB=""
stabv=0
HDF5_ENABLED=""

RED='\033[0;31m'
NC='\033[0m' # No Color

while [ "$#" -gt "0" ]; do
  ARG="$(echo "$1" | tr '[:lower:]' '[:upper:]')"
  shift 1
  case "$ARG" in
    STAB1|STAB2|STAB3|LOG)
      if [ "$stabv" = "1" ]; then
         printf "Additional STAB configuration found. Overwriting %s with %s .\n" "$STAB" "$ARG"
      fi
      STAB="$ARG"
      stabv="1"
    ;;
    NOMPI|MPI|TEMPERING|SERIAL)
      if [ "$modev" = "1" ]; then
         printf "Additional MODE configuration found. Overwriting %s with %s .\n" "$MODE" "$ARG"
      fi
      MODE="$ARG"
      modev="1"
    ;;
    HDF5)
      HDF5_ENABLED="1"
    ;;
    DEVEL|DEVELOPMENT)
      #DEVEL="1"
      GNUOPTFLAGS="$GNUOPTFLAGS $GNUDEVFLAGS"
      INTELOPTFLAGS="$INTELOPTFLAGS $INTELDEVFLAGS"
      PGIOPTFLAGS="$PGIOPTFLAGS $PGIDEVFLAGS"
    ;;
    *)
      if [ "$Machinev" = "1" ]; then
         printf "Additional MACHINE / unrecognized configuration found. Overwriting %s with %s .\n" "$MACHINE" "$ARG"
      fi
      MACHINE="$ARG"
      Machinev="1"
    ;;
  esac
done

printf "\n"

case $MODE in
  NOMPI|SERIAL)
    printf "serial job.\n"
    PROGRAMMCONFIGURATION=""
    INTELCOMPILER="ifort"
    GNUCOMPILER="gfortran"
    MPICOMP=0
  ;;

  TEMPERING)
    printf "Activating parallel tempering.\n"
    printf "This requires also MPI parallization which is set as well.\n"
    PROGRAMMCONFIGURATION="-DMPI -DTEMPERING"
    INTELCOMPILER="mpiifort"
    GNUCOMPILER="mpifort"
    MPICOMP=1
  ;;

  MPI)
    printf "Activating MPI parallization.\n"
    PROGRAMMCONFIGURATION="-DMPI"
    INTELCOMPILER="mpiifort"
    GNUCOMPILER="mpifort"
    MPICOMP=1
  ;;

  *)
    printf "Activating ${RED}MPI parallization (default)${NC}.\n"
    printf "To turn MPI off, pass noMPI as the second argument.\n"
    printf "To turn on parallel tempering, pass Tempering as the second argument.\n"
    PROGRAMMCONFIGURATION="-DMPI"
    INTELCOMPILER="mpiifort"
    GNUCOMPILER="mpifort"
    MPICOMP=1
  ;;
esac

printf "\n"

case $STAB in
  STAB1)
    STABCONFIGURATION="${STABCONFIGURATION} -DSTAB1"
    printf "Using older stabilization with UDV decompositions\n"
  ;;

  STAB2)
    STABCONFIGURATION="${STABCONFIGURATION} -DSTAB2"
    printf "Using older stabilization with UDV decompositions and additional normalizations\n"
  ;;

  STAB3)
    STABCONFIGURATION="${STABCONFIGURATION} -DSTAB3"
    printf "Using newest stabilization which seperates large and small scales\n"
  ;;

  LOG)
    STABCONFIGURATION="${STABCONFIGURATION} -DLOG"
    printf "Using log storage for internal scales\n"
  ;;

  *)
    printf "Using ${RED}default stabilization${NC}\n"
    printf "Possible alternative options are STAB1, STAB2, STAB3 and LOG\n"
  ;;
esac

case $MACHINE in
  #Fakhers MacBook
  FAKHERSMAC)
    # F90OPTFLAGS=$GNUOPTFLAGS
    F90OPTFLAGS="$GNUOPTFLAGS -Wconversion  -Wuninitialized  -fcheck=all -g -fbacktrace"
    F90USEFULFLAGS="$GNUUSEFULFLAGS"
    if [ "$MPICOMP" -eq "0" ]; then
    ALF_FC="gfortran"
    else
    ALF_FC="$mpif90"
    fi
    LIB_BLAS_LAPACK="-llapack -lblas -fopenmp"
    if [ "${HDF5_ENABLED}" = "1" ]; then
      set_hdf5_flags gcc gfortran g++ || return 1
    fi
  ;;

  #LRZ enviroment
  SUPERMUC)
    module switch mpi.ibm  mpi.intel/2018
    module switch intel intel/18.0
    module switch mkl mkl/2018
    module load hdf5

    F90OPTFLAGS="$INTELOPTFLAGS"
    F90USEFULFLAGS="$INTELUSEFULFLAGS"
    ALF_FC="mpiifort"
    LIB_BLAS_LAPACK="$MKL_LIB"
    LIB_HDF5="$HDF5_F90_LIB $HDF5_LIB $SZIP_LIB -lz"
    INC_HDF5="$HDF5_INC"
  ;;

  #LRZ enviroment
  SUPERMUC-NG|NG)
    module switch mpi.intel  mpi.intel/2019
    module switch intel intel/19.0
    module switch mkl mkl/2019
    module load hdf5/serial/1.8
    printf "\n${RED}   !!   unsetting  FORT_BLOCKSIZE  !!${NC}\n"
    unset FORT_BLOCKSIZE

    #module load  mpi.intel
    #module load intel
    #module load mkl

    F90OPTFLAGS="$INTELOPTFLAGS"
    F90USEFULFLAGS="$INTELUSEFULFLAGS"
    ALF_FC="mpiifort"
    LIB_BLAS_LAPACK="$MKL_LIB"
    LIB_HDF5="$HDF5_F90_LIB $HDF5_LIB $SZIP_LIB -lz"
    INC_HDF5="$HDF5_INC"
  ;;

  #JUWELS enviroment
  JUWELS)
    module load Intel
    module load IntelMPI
    module load imkl

    F90OPTFLAGS="$INTELOPTFLAGS"
    F90USEFULFLAGS="$INTELUSEFULFLAGS"
    ALF_FC="mpiifort"
    LIB_BLAS_LAPACK="-mkl"
  ;;

  #Intel (as Hybrid code)
  INTEL)
    F90OPTFLAGS="$INTELOPTFLAGS"
    F90USEFULFLAGS="$INTELUSEFULFLAGS"
    ALF_FC="$INTELCOMPILER"
    LIB_BLAS_LAPACK="-mkl"
    if [ "${HDF5_ENABLED}" = "1" ]; then
      set_hdf5_flags icc ifort icpc || return 1
    fi
  ;;

  #GNU (as Hybrid code)
  GNU)
    F90OPTFLAGS="$GNUOPTFLAGS"
    F90USEFULFLAGS="$GNUUSEFULFLAGS"
    ALF_FC="$GNUCOMPILER"
    LIB_BLAS_LAPACK="-llapack -lblas -fopenmp"
    if [ "${HDF5_ENABLED}" = "1" ]; then
      set_hdf5_flags gcc gfortran g++ || return 1
    fi
  ;;

  #PGI
  PGI)
    F90OPTFLAGS="$PGIOPTFLAGS"
    F90USEFULFLAGS="$PGIUSEFULFLAGS"
    if [ "$MPICOMP" -eq "0" ]; then
      ALF_FC="pgfortran"
    else
      ALF_FC="mpifort"
      printf "\n${RED}   !! Compiler set to 'mpifort' !!\n"
      printf "If this is not your PGI MPI compiler you have to set it manually through:\n"
      printf "    'export ALF_FC=<mpicompiler>'${NC}\n"
    fi
    LIB_BLAS_LAPACK="-llapack -lblas"
    if [ "${HDF5_ENABLED}" = "1" ]; then
      set_hdf5_flags pgcc pgfortran pgc++ || return 1
    fi
  ;;

  #Default (unknown machine)
  *)
    printf "\n"
    printf "${RED}   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!${NC}\n"
    printf "${RED}   !!               UNKNOW MACHINE               !!${NC}\n"
    printf "${RED}   !!         IGNORING PARALLEL SETTINGS         !!${NC}\n"
    printf "${RED}   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!${NC}\n"
    printf "\n"
    printf "Activating fallback option with gfortran for SERIAL JOB.\n"
    printf "\n"
    printf "usage 'source configureHPC.sh MACHINE MODE STAB'\n"
    printf "\n"
    printf "Please choose one of the following machines:\n"
    printf " * SuperMUC\n"
    printf " * SuperMUC-NG\n"
    printf " * JUWELS\n"
    printf " * Devel\n"
    printf " * Intel\n"
    printf " * GNU\n"
    printf " * FakhersMAC\n"
    printf "Possible modes are MPI (default), noMPI and Tempering\n"
    printf "Possible stab are no-argument (default), STAB1 (old), STAB2 (old), STAB3 (newest)\n"
    printf "and LOG (increases accessible scales, e.g. in beta or interaction strength by solving NaN issues)\n"
    printf "Further options: Devel and HDF5"
    printf "To hand an additional flag to the compiler, export it in the varible ALF_FLAGS_EXT prior to soucing this script."

    PROGRAMMCONFIGURATION=""
    F90OPTFLAGS="-cpp -O3 -ffree-line-length-none -ffast-math"
    F90USEFULFLAGS=""

    ALF_FC="gfortran"
    LIB_BLAS_LAPACK="-llapack -lblas"
    if [ "${HDF5_ENABLED}" = "1" ]; then
      set_hdf5_flags gcc gfortran g++ || return 1
    fi
  ;;
esac

PROGRAMMCONFIGURATION="$STABCONFIGURATION $PROGRAMMCONFIGURATION"

Libs="$ALF_DIR/Libraries"
ALF_INC="-I${Libs}/Modules"
ALF_LIB="${Libs}/Modules/modules_90.a ${LIB_BLAS_LAPACK} ${Libs}/libqrref/libqrref.a"
if [ "${HDF5_ENABLED}" = "1" ]; then
  echo; echo "HDF5 enabled"
  ALF_INC="${ALF_INC} ${INC_HDF5}"
  ALF_LIB="${ALF_LIB} ${LIB_HDF5}"
else
  echo; echo "HDF5 disabled"
fi
export ALF_LIB

export ALF_DIR
export ALF_FC

if [ ! -z "${ALF_FLAGS_EXT+x}" ]; then
  printf "\nAppending additional compiler flag '%s'\n" "${ALF_FLAGS_EXT}"
fi

ALF_FLAGS_QRREF="${F90OPTFLAGS} ${ALF_FLAGS_EXT}"
#Modules need to know the programm configuration since entanglement needs MPI
ALF_FLAGS_MODULES="${F90OPTFLAGS} ${PROGRAMMCONFIGURATION} ${ALF_FLAGS_EXT}"
ALF_FLAGS_ANA="${F90USEFULFLAGS} ${F90OPTFLAGS} ${ALF_INC} ${ALF_FLAGS_EXT}"
ALF_FLAGS_PROG="${F90USEFULFLAGS} ${F90OPTFLAGS} ${PROGRAMMCONFIGURATION} ${ALF_INC} ${ALF_FLAGS_EXT}"
# Control with flags -DHDF5 -DHDF5_ZLIB -DOBS_ASCII -DOBS_LEGACY, which observable format to use
if [ "${HDF5_ENABLED}" = "1" ]; then
  ALF_FLAGS_MODULES="${ALF_FLAGS_MODULES} ${INC_HDF5} -DHDF5 -DHDF5_ZLIB"
  ALF_FLAGS_ANA="${ALF_FLAGS_ANA} ${INC_HDF5} -DHDF5 -DHDF5_ZLIB"
  ALF_FLAGS_PROG="${ALF_FLAGS_PROG} -DHDF5 -DHDF5_ZLIB"
fi
export ALF_FLAGS_QRREF
export ALF_FLAGS_MODULES
export ALF_FLAGS_ANA
export ALF_FLAGS_PROG

printf "\nTo compile your program use:    'make TARGET'\n\n"
