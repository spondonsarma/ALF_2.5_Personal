# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/stafusa/ALF/ALF/testsuite

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/stafusa/ALF/ALF/testsuite/tests

# Include any dependencies generated for this target.
include Prog.tests/CMakeFiles/13-Op-Wrapup.dir/depend.make

# Include the progress variables for this target.
include Prog.tests/CMakeFiles/13-Op-Wrapup.dir/progress.make

# Include the compile flags for this target's objects.
include Prog.tests/CMakeFiles/13-Op-Wrapup.dir/flags.make

Prog.tests/CMakeFiles/13-Op-Wrapup.dir/13-Op-Wrapup.F90.o: Prog.tests/CMakeFiles/13-Op-Wrapup.dir/flags.make
Prog.tests/CMakeFiles/13-Op-Wrapup.dir/13-Op-Wrapup.F90.o: ../Prog.tests/13-Op-Wrapup.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/stafusa/ALF/ALF/testsuite/tests/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object Prog.tests/CMakeFiles/13-Op-Wrapup.dir/13-Op-Wrapup.F90.o"
	cd /home/stafusa/ALF/ALF/testsuite/tests/Prog.tests && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/stafusa/ALF/ALF/testsuite/Prog.tests/13-Op-Wrapup.F90 -o CMakeFiles/13-Op-Wrapup.dir/13-Op-Wrapup.F90.o

Prog.tests/CMakeFiles/13-Op-Wrapup.dir/13-Op-Wrapup.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/13-Op-Wrapup.dir/13-Op-Wrapup.F90.i"
	cd /home/stafusa/ALF/ALF/testsuite/tests/Prog.tests && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/stafusa/ALF/ALF/testsuite/Prog.tests/13-Op-Wrapup.F90 > CMakeFiles/13-Op-Wrapup.dir/13-Op-Wrapup.F90.i

Prog.tests/CMakeFiles/13-Op-Wrapup.dir/13-Op-Wrapup.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/13-Op-Wrapup.dir/13-Op-Wrapup.F90.s"
	cd /home/stafusa/ALF/ALF/testsuite/tests/Prog.tests && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/stafusa/ALF/ALF/testsuite/Prog.tests/13-Op-Wrapup.F90 -o CMakeFiles/13-Op-Wrapup.dir/13-Op-Wrapup.F90.s

# Object files for target 13-Op-Wrapup
13__Op__Wrapup_OBJECTS = \
"CMakeFiles/13-Op-Wrapup.dir/13-Op-Wrapup.F90.o"

# External object files for target 13-Op-Wrapup
13__Op__Wrapup_EXTERNAL_OBJECTS =

Prog.tests/13-Op-Wrapup: Prog.tests/CMakeFiles/13-Op-Wrapup.dir/13-Op-Wrapup.F90.o
Prog.tests/13-Op-Wrapup: Prog.tests/CMakeFiles/13-Op-Wrapup.dir/build.make
Prog.tests/13-Op-Wrapup: ../Prog.tests/../../Prog/Operator_mod.o
Prog.tests/13-Op-Wrapup: ../Prog.tests/../../Prog/Fields_mod.o
Prog.tests/13-Op-Wrapup: ../Prog.tests/../../Libraries/Modules/modules_90.a
Prog.tests/13-Op-Wrapup: ../Prog.tests/../../Libraries/libqrref/libqrref.a
Prog.tests/13-Op-Wrapup: /usr/lib/x86_64-linux-gnu/libopenblas.so
Prog.tests/13-Op-Wrapup: /opt/intel/compilers_and_libraries_2018.3.222/linux/mkl/lib/intel64_lin/libmkl_intel_lp64.so
Prog.tests/13-Op-Wrapup: /opt/intel/compilers_and_libraries_2018.3.222/linux/mkl/lib/intel64_lin/libmkl_intel_thread.so
Prog.tests/13-Op-Wrapup: /opt/intel/compilers_and_libraries_2018.3.222/linux/mkl/lib/intel64_lin/libmkl_core.so
Prog.tests/13-Op-Wrapup: /opt/intel/compilers_and_libraries_2018.3.222/linux/compiler/lib/intel64/libiomp5.so
Prog.tests/13-Op-Wrapup: Prog.tests/CMakeFiles/13-Op-Wrapup.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/stafusa/ALF/ALF/testsuite/tests/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable 13-Op-Wrapup"
	cd /home/stafusa/ALF/ALF/testsuite/tests/Prog.tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/13-Op-Wrapup.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Prog.tests/CMakeFiles/13-Op-Wrapup.dir/build: Prog.tests/13-Op-Wrapup

.PHONY : Prog.tests/CMakeFiles/13-Op-Wrapup.dir/build

Prog.tests/CMakeFiles/13-Op-Wrapup.dir/clean:
	cd /home/stafusa/ALF/ALF/testsuite/tests/Prog.tests && $(CMAKE_COMMAND) -P CMakeFiles/13-Op-Wrapup.dir/cmake_clean.cmake
.PHONY : Prog.tests/CMakeFiles/13-Op-Wrapup.dir/clean

Prog.tests/CMakeFiles/13-Op-Wrapup.dir/depend:
	cd /home/stafusa/ALF/ALF/testsuite/tests && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/stafusa/ALF/ALF/testsuite /home/stafusa/ALF/ALF/testsuite/Prog.tests /home/stafusa/ALF/ALF/testsuite/tests /home/stafusa/ALF/ALF/testsuite/tests/Prog.tests /home/stafusa/ALF/ALF/testsuite/tests/Prog.tests/CMakeFiles/13-Op-Wrapup.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Prog.tests/CMakeFiles/13-Op-Wrapup.dir/depend

