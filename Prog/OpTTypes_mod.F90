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


module OpTTypes_mod
    use ContainerElementBase_mod
    implicit none

    type, extends(ContainerElementBase) :: RealOpT
        Real(kind=kind(0.d0)), allocatable, dimension(:,:) :: mat
    contains
        procedure :: simt => RealOpT_simt
        procedure :: mult => RealOpT_mult
    end type RealOpT

    type, extends(ContainerElementBase) :: CmplxOpT
        Complex(kind=kind(0.d0)),allocatable, dimension(:,:):: mat
    contains
        procedure :: simt => CmplxOpT_simt
        procedure :: mult => CmplxOpT_mult
    end type CmplxOpT

contains
    subroutine RealOpT_simt(this, arg)
        class(RealOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), allocatable, dimension(:,:) :: arg
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: temp

    end subroutine
    
    subroutine RealOpT_mult(this, arg)
        class(RealOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), allocatable, dimension(:,:) :: arg
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: temp

    end subroutine

    subroutine CmplxOpT_simt(this, arg)
        class(CmplxOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), allocatable, dimension(:,:) :: arg
    end subroutine
    
    subroutine CmplxOpT_mult(this, arg)
        class(CmplxOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), allocatable, dimension(:,:) :: arg
    end subroutine

end module OpTTypes_mod
