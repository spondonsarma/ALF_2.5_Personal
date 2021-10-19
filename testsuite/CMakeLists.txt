cmake_minimum_required(VERSION 2.8.11)
Project(matmod.tests C Fortran)

#set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran")

add_subdirectory(matmod.tests)
add_subdirectory(Prog.tests)
enable_testing()
