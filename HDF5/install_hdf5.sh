#!/usr/bin/env bash
# Script for automatically downloading and installing HDF5 in current directory
# Needs the following environment variables:
#   CC: C compiler
#   FC: Fortran compiler
#   CXX: C++ compiler
#   HDF5_DIR: Diretory, in which HDF5 gets installed

if [ -d "$HDF5_DIR" ]; then
  printf "\e[31mDirectory %s already exists, aborting HDF5 installation.\e[0m\n" "$HDF5_DIR"
  exit 1
fi

# Create temporary directory
tmpdir=$(mktemp -d 2>/dev/null || mktemp -d -t 'tmpdir')
printf "\e[31mTemporary directory %s created\e[0m\n" "$tmpdir"
cd "$tmpdir" || exit 1

printf "\e[31m========== Downloading source ==========\e[0m\n"
curl https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.7/src/hdf5-1.10.7.tar.gz | tar xz
source_dir="hdf5-1.10.7"

export CC FC CXX
printf "\e[31m=== Build with the following compilers C: %s, Fortran: %s, C++: %s \e[0m\n" "$CC" "$FC" "$CXX"

if ! command -v "$CC" > /dev/null; then
  printf "\e[31m==== C compiler <%s> not available =====\e[0m\n" "$CC"
  exit 1
fi
if ! command -v "$FC" > /dev/null; then
  printf "\e[31m==== FORTRAN compiler <%s> not available =====\e[0m\n" "$FC"
  exit 1
fi
if ! command -v "$CXX" > /dev/null; then
  printf "\e[31m==== C++ compiler <%s> not available =====\e[0m\n" "$CXX"
  exit 1
fi

"$source_dir/configure" --prefix="$HDF5_DIR" --enable-fortran --enable-shared=no --enable-tests=no
if ! make; then
  printf "\e[31m=== Compilation with compilers %s %s in directory %s failed ===\e[0m\n" "$CC" "$FC" "$PWD"
  rm -rf "$HDF5_DIR"
  exit 1
fi
#make check
if ! make install; then
  printf "\e[31m=== Installation of HDF5 in directory %s failed ===\e[0m\n" "$HDF5_DIR"
  rm -rf "$HDF5_DIR"
  exit 1
fi
#make check-install

printf "\e[31mYou can delete the temporary directory %s\e[0m\n" "$tmpdir"
