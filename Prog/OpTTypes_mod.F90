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
    use Operator_mod
    implicit none
    
    type, extends(ContainerElementBase) :: RealOpT
        Real(kind=kind(0.d0)), allocatable, dimension(:,:) :: mat, invmat
        Real(kind=kind(0.d0)) :: g, Zero
        integer, pointer :: P(:)
        Integer :: m, n, Ndim_hop
        
    contains
        procedure :: init => RealOpT_init
        procedure :: simt => RealOpT_simt
        procedure :: rmult => RealOpT_rmult
        procedure :: lmult => RealOpT_lmult
        procedure :: rmultinv => RealOpT_rmultinv
        procedure :: lmultinv => RealOpT_lmultinv
        procedure :: dump => RealOpT_dump
    end type RealOpT

    type, extends(ContainerElementBase) :: CmplxOpT
        Complex(kind=kind(0.d0)),allocatable, dimension(:,:):: mat, invmat
        Complex(kind=kind(0.d0)) :: g
        Real(kind=kind(0.d0)) :: Zero
        integer, pointer :: P(:)
        Integer :: m, n, Ndim_hop
    contains
        procedure :: init => CmplxOpT_init
        procedure :: simt => CmplxOpT_simt
        procedure :: rmult => CmplxOpT_rmult
        procedure :: lmult => CmplxOpT_lmult
        procedure :: rmultinv => CmplxOpT_rmultinv
        procedure :: lmultinv => CmplxOpT_lmultinv
        procedure :: dump => CmplxOpT_dump
    end type CmplxOpT

contains
    subroutine RealOpT_init(this, OpT)
        class(RealOpT) :: this
        Type(Operator), intent(in) :: OpT
        Integer :: i,j
        
        this%Zero = 1.E-12
    end subroutine
    
    subroutine RealOpT_simt(this, arg)
        class(RealOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), allocatable, dimension(:,:) :: arg
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: temp

    end subroutine
    
    subroutine RealOpT_rmult(this, arg)
        class(RealOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout),  dimension(:,:) :: arg
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: out
        Real(kind=kind(0.D0)), allocatable, dimension(:) :: rwork
        Integer :: i, j, sz1, sz2
        
        sz1 = size(arg, 1)
        sz2 = size(arg, 2)
        allocate(out(sz1, sz2), rwork(2*sz1*sz2))
        call zlacrm(sz1, sz2, arg, sz1, this%mat, this%m, out, this%m, rwork) ! zlarcm assumes mat to be square
        arg = out
        deallocate(out, rwork)
    end subroutine
    
    subroutine RealOpT_rmultinv(this, arg)
        class(RealOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout),  dimension(:,:) :: arg
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: out
        Real(kind=kind(0.D0)), allocatable, dimension(:) :: rwork
        Integer :: i, j, sz1, sz2
        
        sz1 = size(arg, 1)
        sz2 = size(arg, 2)
        allocate(out(sz1, sz2), rwork(2*sz1*sz2))
        call zlacrm(sz1, sz2, arg, sz1, this%mat, this%m, out, this%m, rwork) ! zlarcm assumes mat to be square
        arg = out
        deallocate(out, rwork)
    end subroutine
    
    subroutine RealOpT_lmult(this, arg)
        class(RealOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout),  dimension(:,:) :: arg
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: temp

    end subroutine
    
    subroutine RealOpT_lmultinv(this, arg)
        class(RealOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: temp

    end subroutine
    
    subroutine CmplxOpT_init(this, Op_T)
        class(CmplxOpT) :: this
        Type(Operator), intent(in) :: Op_T
        Integer :: i, j
        
        this%Zero = 1.E-12
        this%Ndim_hop = Op_T%N
        this%g = -Op_T%g
!        if (allocated(this%mat) .and. allocated(this%invmat) ) then
            allocate (this%mat(this%Ndim_hop, this%Ndim_hop), this%invmat(this%Ndim_hop, this%Ndim_hop))
!        endif
        Call  Op_exp(this%g, Op_T, this%invmat)
        this%g = Op_T%g
        Call  Op_exp(this%g, Op_T, this%mat )
        DO i = 1, this%Ndim_hop
            DO j = i, this%Ndim_hop
                this%mat(i, j) = (this%mat(i, j) + Conjg(this%mat(j, i)))/2.D0
                this%invmat(i, j) = (this%invmat(i, j) + Conjg(this%invmat(j, i)))/2.D0
            ENDDO
        ENDDO
        this%P => Op_T%P

    end subroutine

    subroutine CmplxOpT_simt(this, arg)
        class(CmplxOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), allocatable, dimension(:,:) :: arg
    end subroutine
    
    subroutine CmplxOpT_rmult(this, arg)
        class(CmplxOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer :: i, j, n1, n2
        
        ! taken from mmthl
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call ZSLHEMM('R', 'U', this%Ndim_hop, n1, n2, this%mat, this%P, arg)
        Endif
    end subroutine
    
        subroutine CmplxOpT_rmultinv(this, arg)
        class(CmplxOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer :: n1, n2
        
        ! taken from mmthl_m1
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call ZSLHEMM('R', 'U', this%Ndim_hop, n1, n2, this%invmat, this%P, arg)
        Endif
    end subroutine
    
    subroutine CmplxOpT_lmult(this, arg)
        class(CmplxOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        integer :: n1, n2
        
        ! taken from mmthr
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call ZSLHEMM('L', 'U', this%Ndim_hop, n1, n2, this%mat, this%P, arg)
        Endif
    end subroutine

    subroutine CmplxOpT_dump(this)
        class(CmplxOpT), intent(in) :: this
        integer :: i,j

        do i = 1, size(this%mat, 1)
write (*,*) (dble(this%mat(i,j)), j = 1,size(this%mat,2) )
enddo
write (*,*) "---------------"
        do i = 1, size(this%mat, 1)
write (*,*) (dble(this%invmat(i,j)), j = 1,size(this%mat,2) )
enddo
        
    end subroutine

    subroutine RealOpT_dump(this)
        class(RealOpT), intent(in) :: this
        integer :: i,j

        do i = 1, size(this%mat, 1)
write (*,*) (dble(this%mat(i,j)), j = 1,size(this%mat,2) )
enddo
write (*,*) "---------------"
        do i = 1, size(this%mat, 1)
write (*,*) (dble(this%invmat(i,j)), j = 1,size(this%mat,2) )
enddo
        
    end subroutine
    
    subroutine CmplxOpT_lmultinv(this, arg)
        class(CmplxOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        integer :: n1, n2
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call ZSLHEMM('L', 'U', this%Ndim_hop, n1, n2, this%invmat, this%P, arg)
        Endif
    end subroutine
    
end module OpTTypes_mod
