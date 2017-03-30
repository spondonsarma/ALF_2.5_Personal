!  Copyright (C) 2017 The ALF project
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
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
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

MODULE UDV_State_mod
    IMPLICIT NONE
    TYPE UDV_State
        COMPLEX (Kind=Kind(0.d0)), allocatable :: U(:, :), V(:, :)
        COMPLEX (Kind=Kind(0.d0)), allocatable :: D(:)
        INTEGER :: ndim
! POINTER and ALLOTABLE give different semantics...
        CONTAINS
            PROCEDURE :: alloc => alloc_UDV_state
            PROCEDURE :: dealloc => dealloc_UDV_state
            PROCEDURE :: reset => reset_UDV_state
            PROCEDURE :: assign => assign_UDV_state
            GENERIC :: ASSIGNMENT(=) => assign
    END TYPE UDV_State

CONTAINS

SUBROUTINE alloc_UDV_state(this, t)
    IMPLICIT NONE
    CLASS(UDV_State), INTENT(INOUT) :: this
    INTEGER :: t
    this%ndim = t
    ALLOCATE(this%U(this%ndim, this%ndim), this%V(this%ndim, this%ndim), this%D(this%ndim))
END SUBROUTINE alloc_UDV_state

SUBROUTINE dealloc_UDV_state(this)
    IMPLICIT NONE
    CLASS(UDV_State), INTENT(INOUT) :: this
    DEALLOCATE(this%U, this%V, this%D)
END SUBROUTINE dealloc_UDV_state

SUBROUTINE reset_UDV_state(this)
    IMPLICIT NONE
    CLASS(UDV_State), INTENT(INOUT) :: this
END SUBROUTINE reset_UDV_state

SUBROUTINE assign_UDV_state(this, rhs)
    IMPLICIT NONE
    CLASS(UDV_State), INTENT(INOUT) :: this
    CLASS(UDV_State), INTENT(IN) :: rhs 
END SUBROUTINE assign_UDV_state

END MODULE UDV_State_mod
