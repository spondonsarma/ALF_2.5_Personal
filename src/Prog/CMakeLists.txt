########################################
# Set up how to compile the source files
########################################

# Add the source files
SET(PROG_src ${SRCPROG}/Hamiltonian_Hub.f90
            #${SRCPROG}/Hamiltonian_SPT.f90
            ${SRCPROG}/Hop_mod.f90
            ${SRCPROG}/Operator.f90
            ${SRCPROG}/UDV_WRAP.f90
            ${SRCPROG}/cgr1.f90
            ${SRCPROG}/cgr2.f90
            ${SRCPROG}/cgr2_1.f90
            ${SRCPROG}/cgr2_2.f90
            ${SRCPROG}/control_mod.f90
            ${SRCPROG}/gperp.f90
            ${SRCPROG}/inconfc.f90
            ${SRCPROG}/main.f90
            ${SRCPROG}/nranf.f90
            ${SRCPROG}/outconfc.f90
            ${SRCPROG}/print_bin_mod.f90
            ${SRCPROG}/tau_m.f90
            ${SRCPROG}/upgrade.f90
            ${SRCPROG}/wrapgrdo.f90
            ${SRCPROG}/wrapgrup.f90
            ${SRCPROG}/wrapul.f90
            ${SRCPROG}/wrapur.f90
)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${PROGEXE} ${PROG_src})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# This links our libraries to the program
TARGET_LINK_LIBRARIES(${PROGEXE} ${MODULES} ${MYEIS} ${MYNAG} ${MYLIN})

# Uncomment if you need to link to BLAS and LAPACK
TARGET_LINK_LIBRARIES(${PROGEXE} ${BLAS_LIBRARIES}
                                ${LAPACK_LIBRARIES}
                                ${CMAKE_THREAD_LIBS_INIT})

# Uncomment if you have parallization
IF(USE_OPENMP)
    SET_TARGET_PROPERTIES(${PROGEXE} PROPERTIES
                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ELSEIF(USE_MPI)
    SET_TARGET_PROPERTIES(${PROGEXE} PROPERTIES
                          COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
                          LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${PROGEXE} ${MPI_Fortran_LIBRARIES})
ENDIF(USE_OPENMP)

#####################################
# Tell how to install this executable
#####################################

IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
INSTALL(TARGETS ${PROGEXE} RUNTIME DESTINATION bin)
