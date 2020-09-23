!  Copyright (C) 2020 The ALF project
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
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
!     
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!     
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       https://alf.physik.uni-wuerzburg.de .
!       
!     - We require the preservation of the above copyright notice and this license in all original files.
!     
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
!       the original source files please visit the homepage https://alf.physik.uni-wuerzburg.de .
! 
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.

! helpful Fortran docs:
! https://materials.prace-ri.eu/400/1/advFortranIntro.pdf
! http://www.chem.helsinki.fi/~manninen/fortran2014/7_Object_oriented_features.pdf

module DynamicMatrixArray_mod
    Use ContainerElementBase_mod
    implicit none

    type :: OpTBasePtrWrapper
        class(ContainerElementBase), allocatable :: dat
    end type

    type :: DynamicMatrixArray
        integer :: avamem ! amount of available space
        integer :: tail ! last index
        Type(OpTbasePtrWrapper), allocatable, dimension(:) :: data
    contains
        procedure :: init => DynamicMatrixArray_init
        procedure :: dealloc => DynamicMatrixArray_dealloc
        procedure :: pushback => DynamicMatrixArray_pushback
        procedure :: at => DynamicMatrixArray_at
        procedure :: back => DynamicMatrixArray_back
        procedure :: length => DynamicMatrixArray_length
        ! FIXME: do we need insert?
    end type DynamicMatrixArray

contains

subroutine DynamicMatrixArray_init(this)
    class(DynamicMatrixArray) :: this
    type(OpTbasePtrWrapper) :: temp
    this%tail = 1
    this%avamem = 4096/(STORAGE_SIZE(temp)/8) ! allocate a page of memory ! Note STORAGE_SIZE: F2008, SIZEOF: GCC Extension
    allocate(this%data(this%avamem))
end subroutine DynamicMatrixArray_init

subroutine DynamicMatrixArray_dealloc(this)
    class(DynamicMatrixArray) :: this
    deallocate(this%data)
end subroutine

subroutine DynamicMatrixArray_pushback(this, itm)
    class(DynamicMatrixArray) :: this
    class(ContainerElementBase), intent(in) :: itm ! Type(...) always has to match exactly, class(...) allows for polymorphism
    type(OpTbasePtrWrapper), allocatable, dimension(:) :: temp
    integer :: i

    if (this%tail == this%avamem) then ! check if this still works the same as for plain ints.
        ! reallocate the memory
        write (*,*) "not enough space -> growing."
        call MOVE_ALLOC(this%data, temp)
        allocate(this%data(2*this%avamem))
        do i = 1, this%avamem
            this%data(i) = temp(i)
        enddo
        deallocate(temp)
        this%avamem = 2*this%avamem
    endif
    this%data(this%tail)%dat = itm
    this%tail = this%tail + 1
end subroutine

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!
!> @brief
!> return a pointer to the object stored at position i
!>
!> @param[in]  i the index
!> @param[out] the content stored at the position i
!
!--------------------------------------------------------------------
function DynamicMatrixArray_at(this, pos) result(itm)
    class(DynamicMatrixArray) :: this
    integer, intent(in) :: pos
    class(ContainerElementBase), allocatable :: itm
    itm = this%data(pos)%dat
end function

subroutine DynamicMatrixArray_back(this, itm)
    class(DynamicMatrixArray) :: this
    class(ContainerElementBase), intent(out), allocatable :: itm
    itm = this%data(this%tail-1)%dat
end subroutine

function DynamicMatrixArray_length(this) result(l)
    class(DynamicMatrixArray) :: this
    integer :: l
    l = this%tail-1
end function

end module DynamicMatrixArray_mod
