
!  Copyright (C) 2016 2017 The ALF project
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
!       http://alf.physik.uni-wuerzburg.de .
!       
!     - We require the preservation of the above copyright notice and this license in all original files.
!     
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
! 
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.
 
     SUBROUTINE WRAPUR(NTAU, NTAU1, UDVR)

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Given    B(NTAU,  1 ) =  UR, DR, VR
!> Returns  B(NTAU1, 1 ) =  UR, DR, VR
!> NOTE:    NTAU1 > NTAU.
!
!-------------------------------------------------------------------

        Use Hop_mod
        Use UDV_State_mod
#if defined(STAB2) || defined(STAB1)         
        Use Hamiltonian
        Use UDV_Wrap_mod
        Implicit None

        ! Arguments
        CLASS(UDV_State), intent(inout) :: udvr(N_FL)
        Integer :: NTAU1, NTAU


        ! Working space.
        Complex (Kind=Kind(0.d0)) :: Z_ONE
        COMPLEX (Kind=Kind(0.d0)) :: V1(Ndim,Ndim), TMP(Ndim,Ndim), TMP1(Ndim,Ndim)
        Integer ::NCON, NT, I, J, n, nf
        Real (Kind=Kind(0.d0)) :: X

        NCON = 0  ! Test for UDV ::::  0: Off,  1: On.
        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        
        Do nf = 1,N_FL
!            CALL INITD(TMP,Z_ONE)
           DO NT = NTAU + 1, NTAU1
              !CALL MMULT(TMP1,Exp_T(:,:,nf) ,TMP)
              Call Hop_mod_mmthr(UDVR(nf)%U,nf)
              Do n = 1,Size(Op_V,1)
!                  X = Phi(nsigma(n,nt),Op_V(n,nf)%type)
                 Call Op_mmultR(UDVR(nf)%U,Op_V(n,nf),nsigma(n,nt),Ndim,'n')
              ENDDO
           ENDDO
!            CALL MMULT(TMP1,TMP, udvr(nf)%U)
           DO J = 1,NDim
              DO I = 1,NDim
                 TMP1(I,J) = UDVR(nf)%U(I,J)*udvr(nf)%D(J)
                 TMP(I,J)  = udvr(nf)%V(I,J)
              ENDDO
           ENDDO
           CALL UDV_WRAP_Pivot(TMP1, udvr(nf)%U, udvr(nf)%D, V1,NCON,Ndim,Ndim)
           CALL MMULT(udvr(nf)%V, V1, TMP)
        ENDDO
#else
        Use Operator_mod, only : Phi
        Implicit None

        ! Arguments
        CLASS(UDV_State), intent(inout) :: udvr(N_FL)
        Integer :: NTAU1, NTAU


        ! Working space.
        Complex (Kind=Kind(0.d0)) :: Z_ONE
!         COMPLEX (Kind=Kind(0.d0)), allocatable, dimension(:, :) :: TMP, TMP1
        Integer :: NT, NCON, n, nf
        Real (Kind=Kind(0.d0)) :: X

        NCON = 0  ! Test for UDV ::::  0: Off,  1: On.
!         Allocate (TMP(Ndim,Ndim), TMP1(Ndim,Ndim))
        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        Do nf = 1,N_FL
!            CALL INITD(TMP,Z_ONE)
           DO NT = NTAU + 1, NTAU1
              !CALL MMULT(TMP1,Exp_T(:,:,nf) ,TMP)
              Call Hop_mod_mmthr(UDVR(nf)%U,nf)
              Do n = 1,Size(Op_V,1)
!                  X = Phi(nsigma(n,nt),Op_V(n,nf)%type)
                 Call Op_mmultR(UDVR(nf)%U,Op_V(n,nf),nsigma(n,nt),Ndim,'n')
              ENDDO
           ENDDO

           CALL UDVR(nf)%left_decompose !(UDVR(nf)%U, TMP1, NCON)
        ENDDO
!         deallocate(TMP, TMP1)

#endif
      END SUBROUTINE WRAPUR
