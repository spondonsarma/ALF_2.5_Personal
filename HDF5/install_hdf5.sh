#!/usr/bin/env bash
#Script for automatically downloading and compiling HDF5

build()
{
	export CC="$2" FC="$3" CXX="$4"
	if ! command -v "$CC" > /dev/null; then
		printf "\e[31m==== C compiler <%s> not available =====\e[0m" "$CC"
		return
	fi
	if ! command -v "$FC" > /dev/null; then
		printf "\e[31m==== FORTRAN compiler <%s> not available =====\e[0m" "$FC"
		return
	fi
	if ! command -v "$CXX" > /dev/null; then
		printf "\e[31m==== C++ compiler <%s> not available =====\e[0m" "$CXX"
		return
	fi

	mkdir "$1" || exit 1
	# shellcheck disable=SC2164
	cd "$1"
	"../$source_dir/configure" --prefix="$dir/$1" --enable-fortran
	if ! make; then
		printf "\e[31m====== Compilation with %s compilers failed =======\e[0m" "$1"
		return
	fi
	#make check
	make install
	#make check-install
	cd ..
}

dir="$PWD"

# Create temporary directory
tmpdir=$(mktemp -d 2>/dev/null || mktemp -d -t 'tmpdir')
printf "\e[31mTemporary directory %s created\e[0m\n" "$tmpdir"
cd "$tmpdir" || exit 1

printf "\e[31m========== Downloading source ==========\e[0m\n"
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.6/src/hdf5-1.10.6.tar.gz

printf "\e[31m========== Unzipping source ==========\e[0m\n"
tar xzf hdf5-1.10.6.tar.gz || exit 1
source_dir="hdf5-1.10.6"

printf "\e[31m========== Build with Intel compilers ==========\e[0m\n"
build intel icc ifort icpc

printf "\e[31m========== Build with GNU compilers ==========\e[0m\n"
build gnu gcc gfortran g++

printf "\e[31m========== Build with PGI compilers ==========\e[0m\n"
build pgi pgcc pgfortran pgc++

printf "\e[31mYou can delete the temporary directory %s\e[0m\n" "$tmpdir"
