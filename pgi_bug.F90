! This is a short test program for a bug encountered in older PGI
! Fortran compilers (19.10 and earlier), which cause a segmentation fault.
! Gets compiled and executed in the Makefile.

Module test_mod
    
Implicit none

  PROCEDURE(proc_base), POINTER, save :: proc

contains
  subroutine ham_set()
    implicit none
    proc => proc_base
  end subroutine ham_set

  Subroutine proc_base(array)
    Implicit none
    Real (Kind=Kind(0.d0)), Intent(inout), allocatable :: array(:,:)
    
    ! print*, "trying to access array"
    array(:,:)  = 0.d0
    ! print*, "done"
  end Subroutine proc_base

end Module test_mod

Program test
  Use test_mod
  implicit none
  
  Real(Kind=Kind(0.d0)), allocatable :: array(:,:)
  
  call ham_set()  ! This works not for older pgfortran versions
  !proc => proc_base  ! This works
  
  Allocate(array(5, 5))

  Call proc(array)
end program
