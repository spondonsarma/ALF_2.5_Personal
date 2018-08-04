!  Copyright (C) 2016 - 2018 The ALF project
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
!       http://alf.physik.uni-wuerzburg.de .
!       
!     - We require the preservation of the above copyright notice and this license in all original files.
!     
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
! 
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.
      
      SUBROUTINE WRAPUL(NTAU1, NTAU, UDVL)

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Given    B(LTROT,NTAU1,Nf  ) =  VL, DL, UL
!> Returns  B(LTROT,NTAU, Nf  ) =  VL, DL, UL
!
!--------------------------------------------------------------------

        !NOTE:    NTAU1 > NTAU.
        Use UDV_State_mod
#if defined(STAB2) ||  defined(STAB1) 
        Use Hamiltonian
        Use Hop_mod
        Use UDV_Wrap_mod

        Implicit none

        ! Arguments
        CLASS(UDV_State), intent(inout), allocatable :: UDVL(N_FL)
        Integer :: NTAU1, NTAU


        ! Working space.
        COMPLEX (Kind=Kind(0.d0)) ::  U1(Ndim,Ndim), V1(Ndim,Ndim), TMP(Ndim,Ndim), TMP1(Ndim,Ndim)
        COMPLEX (Kind=Kind(0.d0)) ::  Z_ONE, beta
        Integer :: NT, NCON, n, nf
        Real    (Kind=Kind(0.d0)) ::  X
 

        NCON = 0  ! Test for UDV ::::  0: Off,  1: On.

        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        beta = 0.D0
        Do nf = 1, N_FL
           CALL INITD(TMP,Z_ONE)
           DO NT = NTAU1, NTAU+1 , -1
              Do n = Size(Op_V,1),1,-1
                 Call Op_mmultL(Tmp,Op_V(n,nf),nsigma%f(n,nt),'n')
              enddo
              !CALL MMULT( TMP1,Tmp,Exp_T(:,:,nf) )
              Call  Hop_mod_mmthl (Tmp,nf)
              ! Tmp = Tmp1
           ENDDO
           
           !Carry out U,D,V decomposition.
           CALL ZGEMM('C', 'N', Ndim, UDVL(nf)%N_part, Ndim, Z_ONE, TMP, Ndim, udvl(nf)%U(1, 1), Ndim, beta, TMP1, Ndim)
           if( ALLOCATED(UDVL(nf)%V) ) then
              DO n = 1,UDVL(nf)%N_part
                  TMP1(:, n) = TMP1(:, n) * udvl(nf)%D(n)
              ENDDO
              CALL UDV_WRAP_Pivot(TMP1,udvl(nf)%U,udvl(nf)%D,V1,NCON,Ndim,Ndim)
              CALL ZGEMM('N', 'C', Ndim, Ndim, Ndim, Z_ONE, udvl(nf)%V(1,1), Ndim, V1, Ndim, beta, TMP1, Ndim)
              udvl(nf)%V = TMP1
           else
              CALL UDV_WRAP_Pivot(TMP1(:,1:UDVL(nf)%N_part),udvl(nf)%U,udvl(nf)%D,V1,NCON,Ndim,UDVL(nf)%N_part)
           endif
        ENDDO

#else
        Use Hop_mod
        Implicit none

        ! Arguments
        CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: UDVL
        Integer, intent(in) :: NTAU1, NTAU


        ! Working space.
!         TYPE(UDV_State) :: udvlocal
!         COMPLEX (Kind=Kind(0.d0)), allocatable, dimension(:, :) :: TMP, TMP1
!         COMPLEX (Kind=Kind(0.d0)) ::  Z_ONE
        Integer :: NT, n, nf
!         Real    (Kind=Kind(0.d0)) ::  X, XMAX, XMEAN
 
        Do nf = 1, N_FL
           DO NT = NTAU1, NTAU+1 , -1
              Do n = Size(Op_V,1),1,-1
                 Call Op_mmultR(udvl(nf)%U,Op_V(n,nf),nsigma%f(n,nt),'c')
              enddo
              Call  Hop_mod_mmthlc (udvl(nf)%U,nf)
           ENDDO
           
           !Carry out U,D,V decomposition.
           CALL UDVL(nf)%decompose
        Enddo
#endif
      END SUBROUTINE WRAPUL
      
