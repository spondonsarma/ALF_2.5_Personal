!  Copyright (C) 2016 - 2018 The ALF project
!
!  This file is part of the ALF project.
!
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with ALF. If not, see http://www.gnu.org/licenses/.
!
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de .
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.


!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> This module provides an interface to allocate memory that is shared
!> between different MPI jobs from a single communicator on the same node
!
!--------------------------------------------------------------------

module mpi_shared_memory

    USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER, C_SIZEOF
#ifdef MPI
    use mpi
#endif
    use iso_fortran_env, only: output_unit, error_unit
    Implicit none
  
    integer, private, save :: nodecomm,noderank
    logical, private, save :: initialized=.false.
    
    INTERFACE allocate_shared_memory
      MODULE PROCEDURE allocate_shared_memory_1Dreal,  allocate_shared_memory_2Dreal,  allocate_shared_memory_3Dreal,  allocate_shared_memory_4Dreal, &
                    &  allocate_shared_memory_1Dcmplx, allocate_shared_memory_2Dcmplx, allocate_shared_memory_3Dcmplx, allocate_shared_memory_4Dcmplx
    END INTERFACE
    
    Contains

      subroutine mpi_shared_memory_init(mpi_communiator)
        Implicit none
        integer, intent(in) :: mpi_communiator
        
        integer :: ierr

#ifdef MPI
        CALL MPI_Comm_split_type(mpi_communiator, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, nodecomm,ierr)
        CALL MPI_Comm_rank(nodecomm, noderank,ierr)
        initialized=.true.
#else
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        error stop 1
#endif
  
      end subroutine mpi_shared_memory_init
      
      subroutine allocate_shared_memory_1Dreal(fortran_array, mpi_win, myrank, arrayshape)
        Implicit none
        real (Kind=Kind(0.d0)), POINTER, intent(inout) :: fortran_array(:)
        integer, intent(out) :: mpi_win
        integer, intent(out) :: myrank
        integer, intent(in)  :: arrayshape(1)
        
#ifdef MPI
        INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
        INTEGER :: disp_unit, ierr
        real(Kind=Kind(0.d0)) :: dummy_real_dp
        TYPE(C_PTR) :: baseptr
        
        if (.not. initialized) then
            WRITE(error_unit,*) 'Please initialize the mpi_shared_memory module before allocating the first array'
            error stop 1
        endif
        myrank = noderank
        if (noderank == 0) then
            windowsize = int(arrayshape(1),MPI_ADDRESS_KIND)*int(C_SIZEOF(dummy_real_dp),MPI_ADDRESS_KIND) !*8 for double ! Put the actual data size here
        else
            windowsize = 0_MPI_ADDRESS_KIND
        end if
        disp_unit = 1
        CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL, nodecomm, baseptr, mpi_win, ierr)

        ! Obtain the location of the memory segment
        if (noderank /= 0) then
            CALL MPI_Win_shared_query(mpi_win, 0, windowsize, disp_unit, baseptr, ierr)
        end if

        ! baseptr can now be associated with a Fortran pointer
        ! and thus used to access the shared data
        CALL C_F_POINTER(baseptr, fortran_array,arrayshape)
#else
        mpi_win=-1
        myrank=-1
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        error stop 1
#endif
  
      end subroutine allocate_shared_memory_1Dreal
      
      subroutine allocate_shared_memory_2Dreal(fortran_array, mpi_win, myrank, arrayshape)
        Implicit none
        real (Kind=Kind(0.d0)), POINTER, intent(inout) :: fortran_array(:,:)
        integer, intent(out) :: mpi_win
        integer, intent(out) :: myrank
        integer, intent(in)  :: arrayshape(2)
        
#ifdef MPI
        INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
        INTEGER :: disp_unit, ierr
        real(Kind=Kind(0.d0)) :: dummy_real_dp
        TYPE(C_PTR) :: baseptr
        
        if (.not. initialized) then
            WRITE(error_unit,*) 'Please initialize the mpi_shared_memory module before allocating the first array'
            error stop 1
        endif
        myrank = noderank
        if (noderank == 0) then
            windowsize = int(arrayshape(1)*arrayshape(2),MPI_ADDRESS_KIND)*int(C_SIZEOF(dummy_real_dp),MPI_ADDRESS_KIND)
        else
            windowsize = 0_MPI_ADDRESS_KIND
        end if
        disp_unit = 1
        CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL, nodecomm, baseptr, mpi_win, ierr)

        ! Obtain the location of the memory segment
        if (noderank /= 0) then
            CALL MPI_Win_shared_query(mpi_win, 0, windowsize, disp_unit, baseptr, ierr)
        end if

        ! baseptr can now be associated with a Fortran pointer
        ! and thus used to access the shared data
        CALL C_F_POINTER(baseptr, fortran_array,arrayshape)
#else
        mpi_win=-1
        myrank=-1
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        error stop 1
#endif
  
      end subroutine allocate_shared_memory_2Dreal
      
      subroutine allocate_shared_memory_3Dreal(fortran_array, mpi_win, myrank, arrayshape)
        Implicit none
        real (Kind=Kind(0.d0)), POINTER, intent(inout) :: fortran_array(:,:,:)
        integer, intent(out) :: mpi_win
        integer, intent(out) :: myrank
        integer, intent(in)  :: arrayshape(3)
        
#ifdef MPI
        INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
        INTEGER :: disp_unit, ierr
        real(Kind=Kind(0.d0)) :: dummy_real_dp
        TYPE(C_PTR) :: baseptr
        
        if (.not. initialized) then
            WRITE(error_unit,*) 'Please initialize the mpi_shared_memory module before allocating the first array'
            error stop 1
        endif
        myrank = noderank
        if (noderank == 0) then
            windowsize = int(arrayshape(1)*arrayshape(2)*arrayshape(3),MPI_ADDRESS_KIND)*int(C_SIZEOF(dummy_real_dp),MPI_ADDRESS_KIND)
        else
            windowsize = 0_MPI_ADDRESS_KIND
        end if
        disp_unit = 1
        CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL, nodecomm, baseptr, mpi_win, ierr)

        ! Obtain the location of the memory segment
        if (noderank /= 0) then
            CALL MPI_Win_shared_query(mpi_win, 0, windowsize, disp_unit, baseptr, ierr)
        end if

        ! baseptr can now be associated with a Fortran pointer
        ! and thus used to access the shared data
        CALL C_F_POINTER(baseptr, fortran_array,arrayshape)
#else
        mpi_win=-1
        myrank=-1
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        error stop 1
#endif
  
      end subroutine allocate_shared_memory_3Dreal
      
      subroutine allocate_shared_memory_4Dreal(fortran_array, mpi_win, myrank, arrayshape)
        Implicit none
        real (Kind=Kind(0.d0)), POINTER, intent(inout) :: fortran_array(:,:,:,:)
        integer, intent(out) :: mpi_win
        integer, intent(out) :: myrank
        integer, intent(in)  :: arrayshape(4)
        
#ifdef MPI
        INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
        INTEGER :: disp_unit, ierr, num_elems
        real(Kind=Kind(0.d0)) :: dummy_real_dp
        TYPE(C_PTR) :: baseptr
        
        if (.not. initialized) then
            WRITE(error_unit,*) 'Please initialize the mpi_shared_memory module before allocating the first array'
            error stop 1
        endif
        myrank = noderank
        if (noderank == 0) then
            num_elems=arrayshape(1)*arrayshape(2)*arrayshape(3)*arrayshape(4)
            windowsize = int(num_elems,MPI_ADDRESS_KIND)*int(C_SIZEOF(dummy_real_dp),MPI_ADDRESS_KIND) !*8 for double ! Put the actual data size here
        else
            windowsize = 0_MPI_ADDRESS_KIND
        end if
        disp_unit = 1
        CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL, nodecomm, baseptr, mpi_win, ierr)

        ! Obtain the location of the memory segment
        if (noderank /= 0) then
            CALL MPI_Win_shared_query(mpi_win, 0, windowsize, disp_unit, baseptr, ierr)
        end if

        ! baseptr can now be associated with a Fortran pointer
        ! and thus used to access the shared data
        CALL C_F_POINTER(baseptr, fortran_array,arrayshape)
#else
        mpi_win=-1
        myrank=-1
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        error stop 1
#endif
  
      end subroutine allocate_shared_memory_4Dreal
      
      subroutine allocate_shared_memory_1Dcmplx(fortran_array, mpi_win, myrank, arrayshape)
        Implicit none
        complex (Kind=Kind(0.d0)), POINTER, intent(inout) :: fortran_array(:)
        integer, intent(out) :: mpi_win
        integer, intent(out) :: myrank
        integer, intent(in)  :: arrayshape(1)
        
#ifdef MPI
        INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
        INTEGER :: disp_unit, ierr
        complex(Kind=Kind(0.d0)) :: dummy_real_dp
        TYPE(C_PTR) :: baseptr
        
        if (.not. initialized) then
            WRITE(error_unit,*) 'Please initialize the mpi_shared_memory module before allocating the first array'
            error stop 1
        endif
        myrank = noderank
        if (noderank == 0) then
            windowsize = int(arrayshape(1),MPI_ADDRESS_KIND)*int(C_SIZEOF(dummy_real_dp),MPI_ADDRESS_KIND) !*8 for double ! Put the actual data size here
        else
            windowsize = 0_MPI_ADDRESS_KIND
        end if
        disp_unit = 1
        CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL, nodecomm, baseptr, mpi_win, ierr)

        ! Obtain the location of the memory segment
        if (noderank /= 0) then
            CALL MPI_Win_shared_query(mpi_win, 0, windowsize, disp_unit, baseptr, ierr)
        end if

        ! baseptr can now be associated with a Fortran pointer
        ! and thus used to access the shared data
        CALL C_F_POINTER(baseptr, fortran_array,arrayshape)
#else
        mpi_win=-1
        myrank=-1
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        error stop 1
#endif
  
      end subroutine allocate_shared_memory_1Dcmplx
      
      subroutine allocate_shared_memory_2Dcmplx(fortran_array, mpi_win, myrank, arrayshape)
        Implicit none
        complex (Kind=Kind(0.d0)), POINTER, intent(inout) :: fortran_array(:,:)
        integer, intent(out) :: mpi_win
        integer, intent(out) :: myrank
        integer, intent(in)  :: arrayshape(2)
        
#ifdef MPI
        INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
        INTEGER :: disp_unit, ierr
        complex(Kind=Kind(0.d0)) :: dummy_real_dp
        TYPE(C_PTR) :: baseptr
        
        if (.not. initialized) then
            WRITE(error_unit,*) 'Please initialize the mpi_shared_memory module before allocating the first array'
            error stop 1
        endif
        myrank = noderank
        if (noderank == 0) then
            windowsize = int(arrayshape(1)*arrayshape(2),MPI_ADDRESS_KIND)*int(C_SIZEOF(dummy_real_dp),MPI_ADDRESS_KIND)
        else
            windowsize = 0_MPI_ADDRESS_KIND
        end if
        disp_unit = 1
        CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL, nodecomm, baseptr, mpi_win, ierr)

        ! Obtain the location of the memory segment
        if (noderank /= 0) then
            CALL MPI_Win_shared_query(mpi_win, 0, windowsize, disp_unit, baseptr, ierr)
        end if

        ! baseptr can now be associated with a Fortran pointer
        ! and thus used to access the shared data
        CALL C_F_POINTER(baseptr, fortran_array,arrayshape)
#else
        mpi_win=-1
        myrank=-1
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        error stop 1
#endif
  
      end subroutine allocate_shared_memory_2Dcmplx
      
      subroutine allocate_shared_memory_3Dcmplx(fortran_array, mpi_win, myrank, arrayshape)
        Implicit none
        complex (Kind=Kind(0.d0)), POINTER, intent(inout) :: fortran_array(:,:,:)
        integer, intent(out) :: mpi_win
        integer, intent(out) :: myrank
        integer, intent(in)  :: arrayshape(3)
        
#ifdef MPI
        INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
        INTEGER :: disp_unit, ierr
        complex(Kind=Kind(0.d0)) :: dummy_real_dp
        TYPE(C_PTR) :: baseptr
        
        if (.not. initialized) then
            WRITE(error_unit,*) 'Please initialize the mpi_shared_memory module before allocating the first array'
            error stop 1
        endif
        myrank = noderank
        if (noderank == 0) then
            windowsize = int(arrayshape(1)*arrayshape(2)*arrayshape(3),MPI_ADDRESS_KIND)*int(C_SIZEOF(dummy_real_dp),MPI_ADDRESS_KIND)
        else
            windowsize = 0_MPI_ADDRESS_KIND
        end if
        disp_unit = 1
        CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL, nodecomm, baseptr, mpi_win, ierr)

        ! Obtain the location of the memory segment
        if (noderank /= 0) then
            CALL MPI_Win_shared_query(mpi_win, 0, windowsize, disp_unit, baseptr, ierr)
        end if

        ! baseptr can now be associated with a Fortran pointer
        ! and thus used to access the shared data
        CALL C_F_POINTER(baseptr, fortran_array,arrayshape)
#else
        mpi_win=-1
        myrank=-1
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        error stop 1
#endif
  
      end subroutine allocate_shared_memory_3Dcmplx
      
      subroutine allocate_shared_memory_4Dcmplx(fortran_array, mpi_win, myrank, arrayshape)
        Implicit none
        complex (Kind=Kind(0.d0)), POINTER, intent(inout) :: fortran_array(:,:,:,:)
        integer, intent(out) :: mpi_win
        integer, intent(out) :: myrank
        integer, intent(in)  :: arrayshape(4)
        
#ifdef MPI
        INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
        INTEGER :: disp_unit, ierr, num_elems
        complex(Kind=Kind(0.d0)) :: dummy_real_dp
        TYPE(C_PTR) :: baseptr
        
        if (.not. initialized) then
            WRITE(error_unit,*) 'Please initialize the mpi_shared_memory module before allocating the first array'
            error stop 1
        endif
        myrank = noderank
        if (noderank == 0) then
            num_elems=arrayshape(1)*arrayshape(2)*arrayshape(3)*arrayshape(4)
            windowsize = int(num_elems,MPI_ADDRESS_KIND)*int(C_SIZEOF(dummy_real_dp),MPI_ADDRESS_KIND) !*8 for double ! Put the actual data size here
        else
            windowsize = 0_MPI_ADDRESS_KIND
        end if
        disp_unit = 1
        CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL, nodecomm, baseptr, mpi_win, ierr)

        ! Obtain the location of the memory segment
        if (noderank /= 0) then
            CALL MPI_Win_shared_query(mpi_win, 0, windowsize, disp_unit, baseptr, ierr)
        end if

        ! baseptr can now be associated with a Fortran pointer
        ! and thus used to access the shared data
        CALL C_F_POINTER(baseptr, fortran_array,arrayshape)
#else
        mpi_win=-1
        myrank=-1
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        error stop 1
#endif
  
      end subroutine allocate_shared_memory_4Dcmplx
      
      subroutine deallocate_shared_memory(mpi_win)
        Implicit none
        integer, intent(in) :: mpi_win
        integer :: ierr
        
#ifdef MPI
        external :: MPI_Win_free    ! This seems to be required by gfortran10 with OpenMPI on Fedora33 (should be part of MPI module)
        call MPI_WIN_FENCE(0, mpi_win, ierr)
        call MPI_BARRIER(nodecomm,ierr)
        call MPI_Win_free(mpi_win,ierr)
#else
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        error stop 1
#endif
  
      end subroutine deallocate_shared_memory
      
end module
