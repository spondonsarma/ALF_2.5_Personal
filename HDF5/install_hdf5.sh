#!/usr/bin/env bash
# Script for automatically downloading and installin HDF5 in current directory
# Needs the follwing environment variables:
#   CC: C compiler
#   FC: Fortran compiler
#   ALF_DIR

dir="$PWD"

# Create temporary directory
tmpdir=$(mktemp -d 2>/dev/null || mktemp -d -t 'tmpdir')
printf "\e[31mTemporary directory %s created\e[0m\n" "$tmpdir"
cd "$tmpdir" || exit 1

printf "\e[31m========== Downloading source ==========\e[0m\n"
curl https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.6/src/hdf5-1.10.6.tar.gz | tar xz
source_dir="hdf5-1.10.6"

export CC FC
printf "\e[31m=== Build with the following compilers C: %s, Fortran: %s \e[0m\n" "$CC" "$FC"

if ! command -v "$CC" > /dev/null; then
  printf "\e[31m==== C compiler <%s> not available =====\e[0m\n" "$CC"
  exit 1
fi
if ! command -v "$FC" > /dev/null; then
  printf "\e[31m==== FORTRAN compiler <%s> not available =====\e[0m\n" "$FC"
  exit 1
fi
#if ! command -v "$CXX" > /dev/null; then
#  printf "\e[31m==== C++ compiler <%s> not available =====\e[0m\n" "$CXX"
#  exit 1
#fi

"$source_dir/configure" --prefix="$dir" --enable-fortran
if ! make; then
  printf "\e[31m=== Compilation with compilers %s %s in directory %s failed ===\e[0m\n" "$CC" "$FC" "$PWD"
  rm -r "$dir"
  exit 1
fi
#make check
make install
#make check-install

printf "\e[31mYou can delete the temporary directory %s\e[0m\n" "$tmpdir"
