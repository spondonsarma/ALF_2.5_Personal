!  Copyright (C) 2016 - 2020 The ALF project
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


      Module Langevin_HMC_mod
        
        Use Hamiltonian
        Use UDV_State_mod
        Use Control
        Use Hop_mod

        
        Implicit none
        
        
        Real    (Kind=Kind(0.d0)),  allocatable, private ::  Forces_0(:,:)



      Contains

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!>   Computes the  forces as well as, on demand,  observables.  
!>   On input:  a)  GR is on the first time slice and  the storage is full with left propagations.
!>              b)  Udvl is on time slice 0.
!>   On output.
!>              a)  Forces (only for field type 3 (i.e. continuous fieds)  are computed   and stored in Forces(:,:)
!>              
!> 
!--------------------------------------------------------------------
        
      SUBROUTINE  Langevin_HMC_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst, LOBS_ST, LOBS_EN, Forces )
        Implicit none
        
        Interface
           SUBROUTINE WRAPUR(NTAU, NTAU1, UDVR)
             Use Hamiltonian
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: UDVR
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUR
           SUBROUTINE WRAPUL(NTAU1, NTAU, UDVL)
             Use Hamiltonian
             Use UDV_State_mod
             Implicit none
             CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: UDVL
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUL
           SUBROUTINE CGR(PHASE,NVAR, GRUP, udvr, udvl)
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), INTENT(IN) :: UDVL, UDVR
             COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(Inout) :: GRUP
             COMPLEX(Kind=Kind(0.d0)) :: PHASE
             INTEGER         :: NVAR
           END SUBROUTINE CGR
        end Interface
        
        CLASS(UDV_State), intent(inout), allocatable, dimension(:  ) :: udvl, udvr
        CLASS(UDV_State), intent(in), allocatable, dimension(:,:)    :: udvst
        Complex (Kind=Kind(0.d0)), intent(inout)                     :: Phase
        Complex (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:)   :: Test
        COMPLEX (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:,:) :: GR, GR_Tilde
        Integer, intent(in),  dimension(:), allocatable :: Stab_nt
        Integer, intent(in) :: LOBS_ST, LOBS_EN
        Complex (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:) ::  Forces 
        

        !Local
        Integer :: NSTM, n, nf, NST, NTAU, nt, nt1, Ntau1, NVAR, N_Type, I, J
        Complex (Kind=Kind(0.d0)) :: Z, Z1
        Real    (Kind=Kind(0.d0)) :: spin
        
        NSTM = Size(Stab_nt,1) - 1 
        !Do  n = 0,NSTM
        !   Write(6,*)  n, Stab_nt(n)
        !Enddo
        
        Forces = cmplx(0.d0,0.d0,Kind(0.d0))
        do nf = 1,N_FL
           if (Projector) then
              CALL udvr(nf)%reset('r',WF_R(nf)%P)
           else
              CALL udvr(nf)%reset('r')
           endif
        Enddo
        NST = 1
        DO NTAU = 0, LTROT-1
           NTAU1 = NTAU + 1
           
           Do nf = 1,N_FL
              CALL HOP_MOD_mmthr   (GR(:,:,nf), nf )
              CALL HOP_MOD_mmthl_m1(GR(:,:,nf), nf )
           Enddo
           Do n = 1, size(OP_V,1) 
              Do nf = 1, N_FL
                 spin = nsigma%f(n,ntau1) ! Phi(nsigma(n,ntau1),Op_V(n,nf)%type)
                 N_type = 1
                 Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),spin,Ndim,N_Type)
                 N_type =  2
                 Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),spin,Ndim,N_Type)
              enddo
              Forces(n,ntau1) = Langevin_HMC_Calc_Force(Gr, n )
           enddo

           If (NTAU1 == Stab_nt(NST) ) then 
              NT1 = Stab_nt(NST-1)
              CALL WRAPUR(NT1, NTAU1, udvr)
              Z = cmplx(1.d0, 0.d0, kind(0.D0))
              Do nf = 1, N_FL
                 ! Read from storage left propagation from LTROT to  NTAU1
                 udvl(nf) = udvst(NST, nf)
                 NVAR = 1
                 IF (NTAU1 .GT. LTROT/2) NVAR = 2
                 TEST(:,:) = GR(:,:,nf)
                 CALL CGR(Z1, NVAR, GR(:,:,nf), UDVR(nf), UDVL(nf))
                 Z = Z*Z1
                 Call Control_PrecisionG(GR(:,:,nf),Test,Ndim)
              ENDDO
              call Op_phase(Z,OP_V,Nsigma,N_SUN) 
              Call Control_PrecisionP(Z,Phase)
              Phase = Z
              NST = NST + 1
           ENDIF
           
           IF (NTAU1.GE. LOBS_ST .AND. NTAU1.LE. LOBS_EN ) THEN
              If (Symm) then
                 Call Hop_mod_Symm(GR_Tilde,GR)
                 CALL Obser( GR_Tilde, PHASE, Ntau1 )
              else
                 CALL Obser( GR, PHASE, Ntau1 )
              endif
           endif
        enddo

      end SUBROUTINE Langevin_HMC_Forces

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!>   For a given  Green function  calculated at position n in the sequence of 
!>   operators,  this routine  computes the force.
!>   
!--------------------------------------------------------------------

       Complex (Kind=kind(0.d0)) function  Langevin_HMC_Calc_Force(Gr,n )

        Implicit none
        
        Complex (Kind=Kind(0.d0)), intent(inout), dimension(:,:,:) :: Gr
        Integer                     :: n

        
        !Local
        Complex (Kind=Kind(0.d0)) :: Z, Z1
        Integer ::  nf, I, J
        
        Langevin_HMC_Calc_Force = cmplx(0.d0,0.d0,Kind(0.d0))
        if (OP_V(n,1)%type == 3 ) then
           Do nf = 1, N_Fl
              Z = cmplx(0.d0,0.d0,Kind(0.d0))
              do I = 1,size(OP_V(n,nf)%P,1)
                 do J = 1,size(OP_V(n,nf)%P,1)
                    Z1 =  cmplx(0.d0,0.d0,Kind(0.d0))
                    if ( I == J ) Z1 = cmplx(1.d0,0.d0,Kind(0.d0))
                    Z  = Z +    Op_V(n,nf)%O(I,J) * ( Z1 - Gr(Op_V(n,nf)%P(J),Op_V(n,nf)%P(I), nf) )
                 Enddo
              Enddo
              Langevin_HMC_Calc_Force =  Langevin_HMC_Calc_Force  - &
                   &    Op_V(n,nf)%g * Z *  cmplx(real(N_SUN,Kind(0.d0)), 0.d0, Kind(0.d0)) 
           Enddo
        endif

      end function Langevin_HMC_Calc_Force
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!>   This routine is called after a  Langevin or HMC step.  It fills the storage with left
!>   propagations and calculates the Green function of the zeroth time slice. 
!>   
!--------------------------------------------------------------------
      Subroutine Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)

        Implicit none
        
        Interface
           SUBROUTINE WRAPUR(NTAU, NTAU1, UDVR)
             Use Hamiltonian
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: UDVR
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUR
           SUBROUTINE WRAPUL(NTAU1, NTAU, UDVL)
             Use Hamiltonian
             Use UDV_State_mod
             Implicit none
             CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: UDVL
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUL
           SUBROUTINE CGR(PHASE,NVAR, GRUP, udvr, udvl)
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), INTENT(IN) :: UDVL, UDVR
             COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(Inout) :: GRUP
             COMPLEX(Kind=Kind(0.d0)) :: PHASE
             INTEGER         :: NVAR
           END SUBROUTINE CGR
        end Interface
        
        CLASS(UDV_State), intent(inout), allocatable, dimension(:  ) :: udvl, udvr
        CLASS(UDV_State), intent(inout), allocatable, dimension(:,:) :: udvst
        Complex (Kind=Kind(0.d0)), intent(inout) :: Phase
        COMPLEX (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:,:) :: GR 
        Integer, intent(in),  dimension(:), allocatable :: Stab_nt
                  
        ! Local
        Integer :: NSTM, nf,  nt, nt1,  NST, NVAR
        Complex (Kind=Kind(0.d0)) :: Z

        
        NSTM = Size(Stab_nt,1) - 1 
        Do nf = 1,N_FL
           if (Projector) then
              CALL udvl(nf)%reset('l',WF_L(nf)%P)
              CALL udvst(NSTM, nf)%reset('l',WF_L(nf)%P)
           else
              CALL udvl(nf)%reset('l')
              CALL udvst(NSTM, nf)%reset('l')
           endif
        ENDDO


        DO NST = NSTM-1,1,-1
           NT1 = Stab_nt(NST+1)
           NT  = Stab_nt(NST  )
           !Write(6,*)'Hi', NT1,NT, NST
           CALL WRAPUL(NT1, NT, UDVL)
           Do nf = 1,N_FL
              UDVST(NST, nf) = UDVL(nf)
           ENDDO
        ENDDO
        NT1 = stab_nt(1)
        CALL WRAPUL(NT1, 0, UDVL)
        
        do nf = 1,N_FL
           if (Projector) then
              CALL udvr(nf)%reset('r',WF_R(nf)%P)
           else
              CALL udvr(nf)%reset('r')
           endif
        ENDDO
        
        NVAR = 1
        Phase = cmplx(1.d0, 0.d0, kind(0.D0))
        do nf = 1,N_Fl
           CALL CGR(Z, NVAR, GR(:,:,nf), UDVR(nf), UDVL(nf))
           Phase = Phase*Z
        Enddo
        call Op_phase(Phase,OP_V,Nsigma,N_SUN)

      end Subroutine Langevin_HMC_Reset_storage
      
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Handles a  Langevin sweep.
!>   On input: a) GR is on the first time slice and  the storage is full with
!>                left propagations.   Udvr  and Udvl are on time slice 1.
!>             b) If L_Forces = .T. (.F.) Fermion_Forces  are (not)  provided      
!>                If L_Forces = .F. (.T.) equal time measurements are (not)  carried  out. 
!>   On output: The  field configuration is  updated.  GR, Udvr,  Udvl and Udvst are as on input but with the
!>              updated configuration.  
!> 
!--------------------------------------------------------------------

      SUBROUTINE  Langevin_update(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst, LOBS_ST, LOBS_EN, &
           &       Forces, L_Forces, Delta_t_Langevin_HMC, Max_Force)
        Implicit none
        
        
        CLASS(UDV_State), intent(inout), allocatable, dimension(:  ) :: udvl, udvr
        CLASS(UDV_State), intent(inout), allocatable, dimension(:,:) :: udvst
        Complex (Kind=Kind(0.d0)), intent(inout) :: Phase
        Complex (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:)   :: Test
        COMPLEX (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:,:) :: GR, GR_Tilde
        Integer, intent(in),  dimension(:), allocatable :: Stab_nt
        Integer, intent(in) :: LOBS_ST, LOBS_EN
        Complex (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:) ::  Forces 
        Logical :: L_Forces 
        Real    (Kind=Kind(0.d0)), intent(in) :: Delta_t_Langevin_HMC, Max_Force

        !Local
        Integer :: N_op, n, nt
        Real    (Kind=Kind(0.d0)) :: X, Xmax

        Real    (Kind=Kind(0.d0)) :: Delta_t_running
        

        If ( .not. L_Forces) &
             &  Call Langevin_HMC_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst,&
             &  LOBS_ST, LOBS_EN, Forces )
        
        Call Control_Langevin   ( Forces,Group_Comm )

        Call Ham_Langevin_HMC_S0_Params(Forces_0,Delta_t_running, "Get" )
          
        N_op = size(nsigma%f,1)
        !  Determine running time step
        Xmax = 0.d0
        do n = 1,N_op
           do nt = 1,Ltrot
              X = abs(Real(Forces(n,nt), Kind(0.d0)))
              if (X > Xmax) Xmax = X
           enddo
        enddo
        Delta_t_running = Delta_t_Langevin_HMC 
        If ( Xmax >  Max_Force ) Delta_t_running = Max_Force * Delta_t_Langevin_HMC / Xmax
        
        Call Ham_Langevin_HMC_S0_Params(Forces_0,Delta_t_running,  "Put" )

        do n = 1,N_op
           if (OP_V(n,1)%type == 3 ) then
              do nt = 1,Ltrot
                 nsigma%f(n,nt)   = nsigma%f(n,nt)  -  ( Forces_0(n,nt) +  &
                      &  real( Phase*Forces(n,nt),kind(0.d0)) / Real(Phase,kind(0.d0)) ) * Delta_t_running + &
                      &  sqrt( 2.d0 * Delta_t_running) * rang_wrap()
              enddo
           endif
        enddo
        Call Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)
        L_Forces = .False. 
        
      END SUBROUTINE LANGEVIN_UPDATE
     

      
      SUBROUTINE  Langevin_setup ( Forces )
        Implicit none

        Integer :: Nr,Nt
        Complex (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:) ::  Forces 
        
        Nr = size(nsigma%f,1)
        Nt = size(nsigma%f,2)
        Allocate ( Forces(Nr,Nt),  Forces_0(Nr,Nt) )
        !Write(6,*) "Langevin: ", Nr,Nt
      end SUBROUTINE Langevin_setup


      SUBROUTINE  Langevin_clear ( Forces )
        Implicit none
        Complex (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:) ::  Forces 
        Deallocate ( Forces, Forces_0 )
      end SUBROUTINE Langevin_clear

    end Module Langevin_HMC_mod
