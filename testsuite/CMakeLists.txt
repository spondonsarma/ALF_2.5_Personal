cmake_minimum_required(VERSION 2.8.11)
set(CMAKE_Fortran_COMPILER $ENV{ALF_FC})

Project(matmod.tests C Fortran)

add_subdirectory(matmod.tests)
add_subdirectory(Prog.tests)
enable_testing()
