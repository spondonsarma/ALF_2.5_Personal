#!/bin/sh
# This script sets necessary environment variables for compiling ALF.
# You need to source it prior to executing make.
USAGE="usage 'source configure.sh MACHINE MODE STAB' \n\
    \n\
Please choose one of the following MACHINEs:\n\
 * GNU\n\
 * Intel\n\
 * PGI\n\
 * SuperMUC-NG\n\
 * JUWELS\n\
Possible MODEs are:\n\
 * MPI (default)\n\
 * noMPI\n\
 * Tempering\n\
Possible STABs are:
 * <no-argument> (default)\n\
 * STAB1 (old)\n\
 * STAB2 (old)\n\
 * STAB3 (newest)\n\
 * LOG (increases accessible scales, e.g. in beta or interaction strength by solving NaN issues)\n\
Further optional arguments: \n\
  Devel: Compile with additional flags for development and debugging\n\
  HDF5: Compile with HDF5\n\
  NO-INTERACTIVE: Do not ask for user confirmation during excution of this script\n\
To hand an additional flag to the compiler, export it in the varible ALF_FLAGS_EXT prior to sourcing this script.\n

For more details check the documentation.\n"

STABCONFIGURATION=""
# STABCONFIGURATION="${STABCONFIGURATION} -DQRREF"

export ALF_DIR="$PWD"

set_hdf5_flags()
{
  CC="$1" FC="$2" CXX="$3"
  
  $FC -o get_compiler_version.out get_compiler_version.F90
  compiler_vers=$(./get_compiler_version.out | sed 's/[ ,()]/_/g')
  
  HDF5_DIR="$ALF_DIR/HDF5/$compiler_vers"
  if [ ! -d "$HDF5_DIR" ]; then
    printf "\nHDF5 is not yet installed for this compiler.\n"
    if [ "$NO_INTERACTIVE" = "" ]; then
      printf "Do you want download and install it now locally in the ALF folder? (Y/n):"
      read -r yn
    else
      yn="Y"
    fi
    case "$yn" in
      y|Y|"")
        printf "${RED}Downloading and installing HDF5 in %s.${NC}\n" "$HDF5_DIR"
        CC="$CC" FC="$FC" CXX="$CXX" HDF5_DIR="$HDF5_DIR" "$ALF_DIR/HDF5/install_hdf5.sh" || return 1
      ;;
      *) 
        printf "Skipping installation of HDF5.\n"
        return 1
      ;;
    esac
  fi
  INC_HDF5="-I$HDF5_DIR/include"
  LIB_HDF5="-L$HDF5_DIR/lib $HDF5_DIR/lib/libhdf5hl_fortran.a $HDF5_DIR/lib/libhdf5_hl.a"
  LIB_HDF5="$LIB_HDF5 $HDF5_DIR/lib/libhdf5_fortran.a $HDF5_DIR/lib/libhdf5.a -lz -ldl -lm -Wl,-rpath -Wl,$HDF5_DIR/lib"
}

check_libs()
{
    FC="$1" LIBS="$2"
    if command -v "$FC" > /dev/null; then       # Compiler binary found
        sh -c "$FC check_libs.f90 $LIBS -o check_libs.out"
        if [ $? -eq 0 ]; then                   # Compiling with $LIBS is successful
            ./check_libs.out || (
              printf "${RED}\n==== Error: Execution of test program using compiler <%s> ====${NC}\n" "$FC"
              printf "${RED}==== and linear algebra libraries <%s> not successful. ====${NC}\n\n" "$LIBS"
              return 1
              )
        else
            printf "${RED}\n==== Error: Linear algebra libraries <%s> not found. ====${NC}\n\n" "$LIBS"
            return 1
        fi
    else
        printf "${RED}\n==== Error: Compiler <%s> not found. ====${NC}\n\n" "$FC"
        return 1
    fi
}

check_python()
{
    if ! command -v python3 > /dev/null; then
	printf "${RED}\n==== Error: Python 3 not found. =====${NC}\n\n"
	return 1
    fi
}

# default optimization flags for Intel compiler
INTELOPTFLAGS="-cpp -O3 -fp-model fast=2 -xHost -unroll -finline-functions -ipo -ip -heap-arrays 1024 -no-wrap-margin"
# INTELOPTFLAGS="-cpp -O3 "
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
GNUDEVFLAGS="-Wconversion -Werror -Wno-error=cpp -fcheck=all -g -fbacktrace -fmax-errors=10"
GNUUSEFULFLAGS="-std=f2008"

# default optimization flags for PGI compiler
PGIOPTFLAGS="-Mpreprocess -O1"
# uncomment the next line if you want to use additional openmp parallelization
PGIOPTFLAGS="${PGIOPTFLAGS} -mp"
PGIDEVFLAGS="-Minform=inform -C -g -traceback"
PGIUSEFULFLAGS=""

MACHINE=""
Machinev=0
MODE=""
modev=0
STAB=""
stabv=0
HDF5_ENABLED=""
NO_INTERACTIVE=""

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
    NO-INTERACTIVE)
      NO_INTERACTIVE="1"
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
    STABCONFIGURATION="${STABCONFIGURATION} -DSTABLOG"
    printf "Using log storage for internal scales\n"
  ;;

  *)
    printf "Using ${RED}default stabilization${NC}\n"
    printf "Possible alternative options are STAB1, STAB2, STAB3 and LOG\n"
  ;;
esac

case $MACHINE in
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

  #LRZ enviroment
  SUPERMUC-NG|NG)
    module load hdf5/1.10.7-intel21
    printf "\n${RED}   !!   unsetting  FORT_BLOCKSIZE  !!${NC}\n"
    unset FORT_BLOCKSIZE

    F90OPTFLAGS="$INTELOPTFLAGS"
    F90USEFULFLAGS="$INTELUSEFULFLAGS"
    ALF_FC="mpiifort"
    LIB_BLAS_LAPACK="$MKL_LIB"
    LIB_HDF5="$HDF5_F90_SHLIB $HDF5_SHLIB"
    INC_HDF5="$HDF5_INC"
  ;;

  #JUWELS enviroment
  JUWELS)
    module load Intel
    module load IntelMPI
    module load imkl
    module load HDF5/1.10.6

    F90OPTFLAGS="$INTELOPTFLAGS"
    F90USEFULFLAGS="$INTELUSEFULFLAGS"
    ALF_FC="mpiifort"
    LIB_BLAS_LAPACK="-mkl"
    LIB_HDF5="–lh5df_fortran"
    INC_HDF5=""
  ;;
  #Default (unknown machine)
  *)
    printf "\n"
    printf "${RED}   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!${NC}\n"
    printf "${RED}   !!               UNKNOW MACHINE               !!${NC}\n"
    printf "${RED}   !!         IGNORING PARALLEL SETTINGS         !!${NC}\n"
    printf "${RED}   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!${NC}\n"
    printf "\n"
    printf "Activating fallback option with gfortran for SERIAL JOB - Deactivating MPI.\n"
    printf "\n"
    printf "$USAGE"
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

check_libs "$ALF_FC" "${LIB_BLAS_LAPACK}" || return 1

check_python || return 1

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

if [ -n "${ALF_FLAGS_EXT+x}" ]; then
  printf "\nAppending additional compiler flag '%s'\n" "${ALF_FLAGS_EXT}"
fi

ALF_FLAGS_QRREF="${F90OPTFLAGS} ${ALF_FLAGS_EXT}"
# Modules need to know the programm configuration since entanglement needs MPI
ALF_FLAGS_MODULES="${F90OPTFLAGS} ${PROGRAMMCONFIGURATION} ${ALF_FLAGS_EXT}"
ALF_FLAGS_ANA="${F90USEFULFLAGS} ${F90OPTFLAGS} ${ALF_INC} ${ALF_FLAGS_EXT}"
ALF_FLAGS_PROG="${F90USEFULFLAGS} ${F90OPTFLAGS} ${PROGRAMMCONFIGURATION} ${ALF_INC} ${ALF_FLAGS_EXT}"
# Control with flags -DHDF5 -DHDF5_ZLIB -DOBS_LEGACY, which observable format to use
if [ "${HDF5_ENABLED}" = "1" ]; then
  ALF_FLAGS_MODULES="${ALF_FLAGS_MODULES} ${INC_HDF5} -DHDF5 -DHDF5_ZLIB"
  ALF_FLAGS_ANA="${ALF_FLAGS_ANA} ${INC_HDF5} -DHDF5 -DHDF5_ZLIB"
  ALF_FLAGS_PROG="${ALF_FLAGS_PROG} -DHDF5 -DHDF5_ZLIB"
fi
export ALF_FLAGS_QRREF
export ALF_FLAGS_MODULES
export ALF_FLAGS_ANA
export ALF_FLAGS_PROG

printf "\nTo compile your program use:    'make'\n\n"
