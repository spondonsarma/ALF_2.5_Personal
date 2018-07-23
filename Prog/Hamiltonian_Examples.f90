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
!       http://alf.physik.uni-wuerzburg.de 
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing   
!       to the ALF project or to mark your material in a reasonable way as different from the original version 


!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>
!> @brief 
!> This module defines the  Hamiltonian and observables.  Here, we have included a
!> set of predefined Hamiltonians. They include the Hubbard and SU(N) tV models
!> on honeycomb, pi-flux and square lattices.

!> @details
!> The public variables of this module are the following
!>
!> 
!> @param [public] OP_V
!> \verbatim
!> Type (Operator), dimension(:,:), allocatable 
!> List of operators of type=1,2 and 3 describing the sequence of interactions on a time slice.
!> The first index runs over this sequence. The second corresponds to the flavor index.  \endverbatim
!> 
!> @param [public] OP_T
!> \verbatim
!> Type (Operator), dimension(:,:), allocatable  
!> Sequence of  operators  accounting for the  hopping on a  time slice. This can include  various
!> checkerboard decompositions. The first index runs over this sequence. The second corresponds to
!> the flavor index. \endverbatim
!> *  The progagation reads:
!> \f$ \prod_{\tau} \; \;  \prod_{n=1}^{N_V}e^{V_n(\tau)}  \prod_{n=1}^{N_T}e^{T_n}  \f$.  That is
!> first the hopping and then the potential energy. 
!>
!> @param [public] WF_L   
!> \verbatim Type (WaveFunction), dimension(:),   allocatable
!> Left trial wave function.  \endverbatim
!>
!> @param [public] WF_R
!> \verbatim Type (WaveFunction), dimension(:),   allocatable
!> Right trial wave function.   For both wave functions the index runs over the flavor index. \endverbatim
!>
!> @param [public]  nsigma(:,:) 
!> \verbatim Integer, allocatable
!> Array containing all auxiliary fields. The first index runs through the operator sequence. The second
!> through the time slies.   \endverbatim
!
!> @param [public]  Ndim
!> \verbatim Integer
!> Total number of orbitals. e.g. # unit cells * # orbitals per unit cell.  \endverbatim
!
!> @param [public]  N_FL
!> \verbatim Integer
!> # of flavors.  Propagation is block diagonal in flavors.  \endverbatim
!
!> @param [public]  N_SUN
!> \verbatim Integer
!> # of colors.  Propagation is color independent.  \endverbatim
!> 
!> @param [public] Ltrot
!> \verbatim Integer
!> Available measurment interval in units of Delta Tau. \endverbatim
!>
!> @param [public] Thtrot  
!>  \verbatim Integer
!> Effective projection parameter in units of Delta Tau.  (Only relevant if projective option is turned on) \endverbatim
!>
!> @param [public] Projector
!> \verbatim Logical
!> Flag for projector. If true then the total number of time slices will correspond to Ltrot + 2*Thtrot \endverbatim
!> 
!> @param [public] Group_Comm 
!> \verbatim Integer
!> Defines MPI communicator  \endverbatim
!
!> @param [public] Symm
!> \verbatim Logical  \endverbatim
!> If set to true then the green functions will be symmetrized
!> before being  sent to the Obser, ObserT subroutines. 
!> In particular, the transformation,  \f$ \tilde{G} =  e^{-\Delta \tau T /2 } G e^{\Delta \tau T /2 } \f$
!> will be carried out  and \f$ \tilde{G} \f$  will be sent to the Obser and ObserT subroutines.  Note that
!> if you want to use this  feature, then you have to be sure the hopping and interaction terms are decomposed
!> symmetrically.
!>  
!> You still have to add some docu for the other private variables in this module.      
!>
!--------------------------------------------------------------------

    Module Hamiltonian

      Use Operator_mod
      Use WaveFunction_mod
      Use Lattices_v3 
      Use MyMats 
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      Use Observables
      Use Predefined_structures

      
      Implicit none

     
      Type (Operator),     dimension(:,:), allocatable  :: Op_V 
      Type (Operator),     dimension(:,:), allocatable  :: Op_T
      Type (WaveFunction), dimension(:),   allocatable  :: WF_L
      Type (WaveFunction), dimension(:),   allocatable  :: WF_R
      Integer, allocatable :: nsigma(:,:)
      Integer              :: Ndim
      Integer              :: N_FL
      Integer              :: N_SUN
      Integer              :: Ltrot
      Integer              :: Thtrot 
      Logical              :: Projector
      Integer              :: Group_Comm
      Logical              :: Symm



      Type (Lattice),       private :: Latt 
      Integer,              private :: L1, L2
      real (Kind=Kind(0.d0)),        private :: ham_T , ham_U,  Ham_chem, Ham_h, Ham_J, Ham_xi, XB_X, Phi_X, Ham_tV
      real (Kind=Kind(0.d0)),        private :: XB_Y, Phi_Y
      real (Kind=Kind(0.d0)),        private :: Dtau, Beta, Theta
      Character (len=64),   private :: Model, Lattice_type
      Logical,              private :: One_dimensional, Checkerboard
      Integer,              private :: N_coord, Norb
      Integer, allocatable, private :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell
      



!>    Privat Observables
      Type (Obser_Vec ),  private, dimension(:), allocatable ::   Obs_scal
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_eq
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_tau
      
!>    Storage for the Ising action
      Real (Kind=Kind(0.d0)),  private :: DW_Ising_tau(-1:1), DW_Ising_Space(-1:1)
      Integer,  allocatable ,  private :: L_bond(:,:), L_bond_inv(:,:), Ising_nnlist(:,:)

      
    contains 


      Subroutine Ham_Set
#ifdef MPI
          Use mpi
#endif
          Implicit none

          integer                :: ierr, N_part
          Real (Kind=Kind(0.d0)) :: Degen
          
          ! L1, L2, Lattice_type, List(:,:), Invlist(:,:) -->  Lattice information
          ! Ham_T, Chem, Phi_X, XB_B, Checkerboard, Symm   -->  Hopping
          ! Interaction                              -->  Model
          
          ! Simulation type                          -->  Finite  T or Projection  Symmetrize Trotter. 
          
          NAMELIST /VAR_Lattice/  L1, L2, Lattice_type, Model,  Checkerboard, N_SUN, Phi_X, XB_X, Symm


          
          NAMELIST /VAR_Hubbard/  ham_T, ham_chem, ham_U,  Dtau, Beta, Theta, Projector

          NAMELIST /VAR_Ising/    ham_T, ham_chem, ham_U, Ham_h, Ham_J, Ham_xi, Dtau, Beta, Theta, Projector

          NAMELIST /VAR_t_V/      ham_T, ham_chem, ham_tV, Dtau, Beta, Theta, Projector

#ifdef MPI
          Integer        :: Isize, Irank
#endif
          ! Global "Default" values.
          N_SUN        = 1
          Checkerboard = .false.
          Symm         = .false.
          Phi_X        = 0.d0
          XB_X         = 1.d0
          Phi_Y        = 0.d0
          XB_Y         = 1.d0

#ifdef MPI
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          If (Irank == 0 ) then
#endif
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(*,*) 'unable to open <parameters>',ierr
                STOP
             END IF
             READ(5,NML=VAR_lattice)
             CLOSE(5)
 
#ifdef MPI
          Endif
          CALL MPI_BCAST(L1          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(L2          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(N_SUN       ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Phi_X       ,1  ,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(XB_X        ,1  ,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Checkerboard,1  ,MPI_LOGICAL  , 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Symm        ,1  ,MPI_LOGICAL  , 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
#endif
          
          Call Predefined_Latt(Lattice_type, L1,L2,Norb,N_coord,Ndim, List,Invlist,Latt)

#ifdef MPI
          If (Irank == 0) then
#endif
             Open (Unit = 50,file="info",status="unknown",position="append")
             Write(50,*) '====================================='
             Write(50,*) 'Model is      : ', Model 
             Write(50,*) 'Lattice is    : ', Lattice_type
             Write(50,*) '# of orbitals : ', Ndim
             If (Lattice_type == "Square" .or. Lattice_type == "One_dimensional" ) then
                Write(50,*) 'X_boundary    : ', XB_X
                Write(50,*) 'Flux_X        : ', Phi_X
             Endif
             Write(50,*) 'Checkerboard  : ', Checkerboard
             Write(50,*) 'Symm. decomp  : ', Symm
             Close(50)
#ifdef MPI
          Endif
#endif


          
#ifdef MPI
          If (Irank == 0 ) then
#endif
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             OPEN(Unit = 50,file="info",status="unknown",position="append")
#ifdef MPI
          Endif
#endif

          Projector = .false.
          Theta = 0.d0
          Thtrot = 0
          Select Case (Model)
          Case ("Hubbard_Mz")
             N_FL  = 2
             N_SUN = 1
#ifdef MPI
             If (Irank == 0 ) then
#endif
                READ(5,NML=VAR_Hubbard)
                Ltrot = nint(beta/dtau)
                if (Projector) Thtrot = nint(theta/dtau)
                Ltrot = Ltrot+2*Thtrot
                if (Projector) then
                  Write(50,*) 'Projective version'
                  Write(50,*) 'Theta         : ', Theta
                  Write(50,*) 'Tau_max       : ', beta
                else
                  Write(50,*) 'Finite temperture version'
                  Write(50,*) 'Beta          : ', Beta
                endif
                Write(50,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
                Write(50,*) 'N_SUN         : ', N_SUN
                Write(50,*) 'N_FL          : ', N_FL
                Write(50,*) 't             : ', Ham_T
                Write(50,*) 'Ham_U         : ', Ham_U
                Write(50,*) 'Ham_chem      : ', Ham_chem
#ifdef MPI
             Endif
#endif
#ifdef MPI
             CALL MPI_BCAST(Ltrot    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Thtrot   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Projector,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_chem ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_U    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Dtau     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Beta     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
          Case ("Hubbard_SU2")
             N_FL = 1
             N_SUN = 2
#ifdef MPI
             If (Irank == 0 ) then
#endif
                READ(5,NML=VAR_Hubbard)
                Ltrot = nint(beta/dtau)
                if (Projector) Thtrot = nint(theta/dtau)
                Ltrot = Ltrot+2*Thtrot
                if (Projector) then
                   Write(50,*) 'Projective version'
                   Write(50,*) 'Theta         : ', Theta
                   Write(50,*) 'Tau_max       : ', beta
                else
                   Write(50,*) 'Finite temperture version'
                   Write(50,*) 'Beta          : ', Beta
                endif
                Write(50,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
                Write(50,*) 'N_SUN         : ', N_SUN
                Write(50,*) 'N_FL          : ', N_FL
                Write(50,*) 't             : ', Ham_T
                Write(50,*) 'Ham_U         : ', Ham_U
                Write(50,*) 'Ham_chem      : ', Ham_chem
#ifdef MPI
             Endif
#endif
#ifdef MPI
             CALL MPI_BCAST(Ltrot    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Thtrot   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Projector,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_T    ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_chem ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_U    ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Dtau     ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Beta     ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
#endif
          Case ("Hubbard_SU2_Ising")
             N_FL = 1
             N_SUN = 2
#ifdef MPI
             If (Irank == 0 ) then
#endif
                READ(5,NML=VAR_Ising)
                Ltrot = nint(beta/dtau)
                if (Projector) Thtrot = nint(theta/dtau)
                Ltrot = Ltrot+2*Thtrot
                if (Projector) then
                  Write(50,*) 'Projective version'
                  Write(50,*) 'Theta         : ', Theta
                  Write(50,*) 'Tau_max       : ', beta
                else
                  Write(50,*) 'Finite temperture version'
                  Write(50,*) 'Beta          : ', Beta
                endif
                Write(50,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
                Write(50,*) 'N_SUN         : ', N_SUN
                Write(50,*) 'N_FL          : ', N_FL
                Write(50,*) 't             : ', Ham_T
                Write(50,*) 'Ham_U         : ', Ham_U
                Write(50,*) 'Ham_chem      : ', Ham_chem
                Write(50,*) 'Ham_xi        : ', Ham_xi
                Write(50,*) 'Ham_J         : ', Ham_J
                Write(50,*) 'Ham_h         : ', Ham_h
#ifdef MPI
             Endif
#endif
             If ( Lattice_type == "Honeycomb" .or. Lattice_type == "Pi_Flux" ) then 
                Write(6,*) "Hubbard_SU2_Ising is only implemented for a square lattice"
                Stop
             Endif
#ifdef MPI
             CALL MPI_BCAST(Ltrot    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Thtrot   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Projector,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_chem ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_U    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Dtau     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Beta     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Ham_xi   ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Ham_J    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Ham_h    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
             Call Setup_Ising_action
          Case ("t_V")
             N_FL  = 1
#ifdef MPI
             If (Irank == 0 ) then
#endif
                READ(5,NML=VAR_t_V )
                Ltrot = nint(beta/dtau)
                if (Projector) Thtrot = nint(theta/dtau)
                Ltrot = Ltrot+2*Thtrot
                if (Projector) then
                  Write(50,*) 'Projective version'
                  Write(50,*) 'Theta         : ', Theta
                  Write(50,*) 'Tau_max       : ', beta
                else
                  Write(50,*) 'Finite temperture version'
                  Write(50,*) 'Beta          : ', Beta
                endif
                Write(50,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
                Write(50,*) 'N_SUN         : ', N_SUN
                Write(50,*) 'N_FL          : ', N_FL
                Write(50,*) 't             : ', Ham_T
                Write(50,*) 'Ham_chem      : ', Ham_chem
                Write(50,*) 'Ham_tV         : ', Ham_tV
#ifdef MPI
             Endif
#endif
#ifdef MPI
             CALL MPI_BCAST(Ltrot    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Thtrot   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Projector,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_chem ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_tV   ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Dtau     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Beta     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
          Case default 
             Write(6,*) "Model not yet implemented!"
             Stop
          end Select

          Call Predefined_Hopping(Lattice_type, Norb,N_coord,Ndim, List,Invlist,Latt, &
           &                      Dtau, Ham_T, Ham_Chem, XB_X, XB_Y, Phi_X, Phi_Y, &
           &                      N_FL,  Checkerboard, Symm, OP_T )

          
          
          if (Projector) then
             N_part = Ndim/2
             Call Predefined_TrialWaveFunction(Lattice_type, Norb,N_coord,Ndim,  List,Invlist,Latt, &
                  &                                  N_part, N_FL,  Degen, WF_L, WF_R)

#ifdef MPI
             If (Irank == 0 ) then
#endif
                Write(50,*) 'Degen of trial wave function: ', Degen
#ifdef MPI
             Endif
#endif             
                
          endif
             

! #ifdef MPI
!           If (Irank == 0 )  then
! #endif
! #if defined(STAB1) 
!              Write(50,*) 'STAB1 is defined '
! #endif
! 
! #if defined(QRREF) 
!              Write(50,*) 'QRREF is defined '
! #endif
!              close(50)
!              Close(5)
! #ifdef MPI
!           endif
! #endif

          call Ham_V


        end Subroutine Ham_Set

!===================================================================================           

!===================================================================================           
        Subroutine Ham_V
          
          Implicit none 
          
          Integer :: nf, I, I1, I2,  nc, nc1,  J
          Real (Kind=Kind(0.d0)) :: X

          
          Select case (Model)
          Case ("Hubbard_SU2")  
             !Write(50,*) 'Model is ', Model
             Allocate(Op_V(Ndim,N_FL))
             do nf = 1,N_FL
                do i  = 1, Ndim
                   Call Op_make(Op_V(i,nf),1)
                enddo
             enddo
             Do nf = 1,N_FL
                nc = 0
                Do i = 1,Ndim
                   nc = nc + 1
                   Op_V(nc,nf)%P(1) = I
                   Op_V(nc,nf)%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
                   Op_V(nc,nf)%g      = SQRT(CMPLX(-DTAU*ham_U/(DBLE(N_SUN)), 0.D0, kind(0.D0))) 
                   Op_V(nc,nf)%alpha  = cmplx(-0.5d0,0.d0, kind(0.D0))
                   Op_V(nc,nf)%type   = 2
                   Call Op_set( Op_V(nc,nf) )
                Enddo
             Enddo
          Case ("Hubbard_Mz") 
             Allocate(Op_V(Ndim,N_FL))
             do nf = 1,N_FL
                do i  = 1, Ndim
                   Call Op_make(Op_V(i,nf),1)
                enddo
             enddo
             Do nf = 1,N_FL
                nc = 0
                X = 1.d0
                if (nf == 2) X = -1.d0
                Do i = 1,Ndim
                   nc = nc + 1
                   Op_V(nc,nf)%P(1) = I
                   Op_V(nc,nf)%O(1,1) = cmplx(1.d0, 0.d0, kind(0.D0))
                   Op_V(nc,nf)%g      = X*SQRT(CMPLX(DTAU*ham_U/2.d0, 0.D0, kind(0.D0))) 
                   Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
                   Op_V(nc,nf)%type   = 2
                   Call Op_set( Op_V(nc,nf) )
                Enddo
             Enddo
          Case ("Hubbard_SU2_Ising") 
             Allocate(Op_V(3*Ndim,N_FL))
             do nf = 1,N_FL
                do i  =  1, N_coord*Ndim
                   call Op_make(Op_V(i,nf),2)
                enddo
                do i  = N_coord*Ndim +1 ,  N_coord*Ndim + Ndim ! For Hubbatd
                   Call Op_make(Op_V(i,nf),1)
                enddo
             enddo
             Do nc = 1,Ndim*N_coord   ! Runs over bonds.  Coordination number = 2.
                ! For the square lattice Ndim = Latt%N
                I1 = L_bond_inv(nc,1)
                I2 = I1
                if ( L_bond_inv(nc,2)  == 1 ) I2 = Latt%nnlist(I1,1,0)
                if ( L_bond_inv(nc,2)  == 2 ) I2 = Latt%nnlist(I1,0,1)
                Op_V(nc,1)%P(1) = I1
                Op_V(nc,1)%P(2) = I2
                Op_V(nc,1)%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.D0)) 
                Op_V(nc,1)%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.D0)) 
                Op_V(nc,1)%g = cmplx(-dtau*Ham_xi,0.D0,kind(0.D0))
                Op_V(nc,1)%alpha = cmplx(0d0,0.d0, kind(0.D0)) 
                Op_V(nc,1)%type =1
                Call Op_set( Op_V(nc,1) )
             Enddo
             
             Do i = 1,Ndim
                nc1 = N_coord*Ndim + i
                Op_V(nc1,1)%P(1)   = i
                Op_V(nc1,1)%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
                Op_V(nc1,1)%g      = sqrt(cmplx(-dtau*ham_U/(DBLE(N_SUN)), 0.D0, kind(0.D0)))
                Op_V(nc1,1)%alpha  = cmplx(-0.5d0,0.d0, kind(0.d0))
                Op_V(nc1,1)%type   = 2
                Call Op_set( Op_V(nc1,1) )
             Enddo
          case ("t_V")
             Allocate(Op_V(N_coord*Latt%N,1))
             do i  =  1, N_coord*Latt%N
                call Op_make(Op_V(i,1),2)
             enddo
             select case (Lattice_type)
             case ("Square")
                !Write(6,*) "N_coord, Latt%N", N_coord, Latt%N, Dtau, ham_tV
                nc = 0
                do I = 1,Latt%N
                   I1 = I
                   do nc1 = 1,N_coord
                      nc = nc + 1
                      if (nc1 == 1 ) I2 = latt%nnlist(I,1,0) 
                      if (nc1 == 2 ) I2 = latt%nnlist(I,0,1)
                      Op_V(nc,1)%P(1) = I1
                      Op_V(nc,1)%P(2) = I2
                      Op_V(nc,1)%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.D0)) 
                      Op_V(nc,1)%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
                      Op_V(nc,1)%g      = SQRT(CMPLX(-DTAU*ham_tV, 0.D0, kind(0.D0))) 
                      Op_V(nc,1)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
                      Op_V(nc,1)%type   = 2
                      Call Op_set( Op_V(nc,1) )
                   Enddo
                Enddo
             case ("Pi_Flux")
                nc = 0
                do I = 1,Latt%N
                   I1 = Invlist(I,1)
                   do nc1 = 1,N_coord
                      nc = nc + 1
                      If (nc1 == 1 )  I2 = invlist(I,2)
                      If (nc1 == 2 )  I2 = invlist(Latt%nnlist(I,0, 1),2) 
                      If (nc1 == 3 )  I2 = invlist(Latt%nnlist(I,-1,1),2)
                      If (nc1 == 4 )  I2 = invlist(Latt%nnlist(I,-1,0),2) 
                      Op_V(nc,1)%P(1) = I1
                      Op_V(nc,1)%P(2) = I2
                      Op_V(nc,1)%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.D0)) 
                      Op_V(nc,1)%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
                      Op_V(nc,1)%g     = SQRT(CMPLX(-DTAU*ham_tV, 0.D0, kind(0.D0))) 
                      Op_V(nc,1)%alpha = cmplx(0.d0, 0.d0, kind(0.D0))
                      Op_V(nc,1)%type  = 2
                      Call Op_set( Op_V(nc,1) )
                   Enddo
                Enddo
             case default
                Write(6,*) "Lattice for t-V  not implemented"
                Stop
             End select
          case default
             Write(6,*) "Model not yet implemented!"
             Stop
          end Select
          
        end Subroutine Ham_V

!===================================================================================           
        Real (Kind=Kind(0.d0)) function S0(n,nt)  
          Implicit none
          Integer, Intent(IN) :: n,nt
          Integer :: nt1,I
          S0 = 1.d0
          If ( Op_V(n,1)%type == 1 ) then 
             do i = 1,4
                S0 = S0*DW_Ising_space(nsigma(n,nt)*nsigma(Ising_nnlist(n,i),nt))
             enddo
             nt1 = nt +1 
             if (nt1 > Ltrot) nt1 = 1
             S0 = S0*DW_Ising_tau(nsigma(n,nt)*nsigma(n,nt1))
             nt1 = nt - 1 
             if (nt1 < 1  ) nt1 = Ltrot
             S0 = S0*DW_Ising_tau(nsigma(n,nt)*nsigma(n,nt1))
             If (S0 < 0.d0) Write(6,*) 'S0 : ', S0
          endif
          
        end function S0

!===================================================================================           
        Subroutine Setup_Ising_action
          
          ! This subroutine sets up lists and arrays so as to enable an 
          ! an efficient calculation of  S0(n,nt) 

          Integer :: nc, nth, n, n1, n2, n3, n4, I, I1, n_orientation
          Real (Kind=Kind(0.d0)) :: X_p(2)

          ! Setup list of bonds for the square lattice.
          Allocate (L_Bond(Latt%N,2),  L_bond_inv(Latt%N*N_coord,2) )
          
          nc = 0
          do nth = 1,2*N_coord  
             Do n1= 1, L1/2
                Do n2 = 1,L2
                   nc = nc + 1
                   n_orientation = 1
                   I1 = 1
                   If (nth == 1 ) then
                      X_p = dble(2*n1)*latt%a1_p + dble(n2)*latt%a2_p 
                      I1 = Inv_R(X_p,Latt)
                      n_orientation = 1
                   elseif (nth == 2) then
                      X_p = dble(2*n1)*latt%a1_p + dble(n2)*latt%a2_p  + latt%a1_p
                      I1 = Inv_R(X_p,Latt)
                      n_orientation = 1
                   elseif (nth == 3) then
                      X_p = dble(n2)*latt%a1_p + dble(2*n1)*latt%a2_p 
                      I1 = Inv_R(X_p,Latt)
                      n_orientation = 2
                   elseif (nth == 4) then
                      X_p = dble(n2)*latt%a1_p + dble(2*n1)*latt%a2_p  + latt%a2_p
                      I1 = Inv_R(X_p,Latt)
                      n_orientation = 2
                   endif
                   L_bond(I1,n_orientation) = nc
                   L_bond_inv(nc,1) = I1  
                   L_bond_inv(nc,2) = n_orientation 
                   ! The bond is given by  I1, I1 + a_(n_orientation).
                Enddo
             Enddo
          Enddo
          ! Setup the nearest neigbour lists for the Ising spins. 
          allocate(Ising_nnlist(2*Latt%N,4)) 
          do I  = 1,Latt%N
             n  = L_bond(I,1)
             n1 = L_bond(Latt%nnlist(I, 1, 0),2)
             n2 = L_bond(Latt%nnlist(I, 0, 0),2)
             n3 = L_bond(Latt%nnlist(I, 0,-1),2)
             n4 = L_bond(Latt%nnlist(I, 1,-1),2)
             Ising_nnlist(n,1) = n1
             Ising_nnlist(n,2) = n2
             Ising_nnlist(n,3) = n3
             Ising_nnlist(n,4) = n4
             n  = L_bond(I,2)
             n1 = L_bond(Latt%nnlist(I, 0, 1),1)
             n2 = L_bond(Latt%nnlist(I,-1, 1),1)
             n3 = L_bond(Latt%nnlist(I,-1, 0),1)
             n4 = L_bond(Latt%nnlist(I, 0, 0),1)
             Ising_nnlist(n,1) = n1
             Ising_nnlist(n,2) = n2
             Ising_nnlist(n,3) = n3
             Ising_nnlist(n,4) = n4
          enddo
          DW_Ising_tau  ( 1) = tanh(Dtau*Ham_h)
          DW_Ising_tau  (-1) = 1.D0/DW_Ising_tau(1)
          DW_Ising_Space( 1) = exp(-2.d0*Dtau*Ham_J) 
          DW_Ising_Space(-1) = exp( 2.d0*Dtau*Ham_J) 
          
        End Subroutine Setup_Ising_action
!===================================================================================           
        Subroutine  Alloc_obs(Ltau) 

          Implicit none
          Integer, Intent(In) :: Ltau
          Integer    ::  i, N, Ns,Nt,No
          Character (len=64) ::  Filename

          ! Scalar observables
          Allocate ( Obs_scal(4) )
          Do I = 1,Size(Obs_scal,1)
             select case (I)
             case (1)
                N = 1;   Filename ="Kin"
             case (2)
                N = 1;   Filename ="Pot"
             case (3)
                N = 1;   Filename ="Part"
             case (4)
                N = 1;   Filename ="Ener"
             case default
                Write(6,*) ' Error in Alloc_obs '  
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
          enddo


          ! Equal time correlators
          Allocate ( Obs_eq(5) )
          Do I = 1,Size(Obs_eq,1)
             select case (I)
             case (1)
                Ns = Latt%N;  No = Norb;  Filename ="Green"
             case (2)
                Ns = Latt%N;  No = Norb;  Filename ="SpinZ"
             case (3)
                Ns = Latt%N;  No = Norb;  Filename ="SpinXY"
             case (4)
                Ns = Latt%N;  No = Norb;  Filename ="Den"
             case (5)
                Ns = Latt%N;  No = Norb;  Filename ="SpinT"
             case default
                Write(6,*) ' Error in Alloc_obs '  
             end select
             Nt = 1
             Call Obser_Latt_make(Obs_eq(I),Ns,Nt,No,Filename)
          enddo

          If (Ltau == 1) then 
             ! Equal time correlators
             Allocate ( Obs_tau(4) )
             Do I = 1,Size(Obs_tau,1)
                select case (I)
                case (1)
                   Ns = Latt%N; No = Norb;  Filename ="Green"
                case (2)
                   Ns = Latt%N; No = Norb;  Filename ="SpinZ"
                case (3)
                   Ns = Latt%N; No = Norb;  Filename ="SpinXY"
                case (4)
                   Ns = Latt%N; No = Norb;  Filename ="Den"
                case default
                   Write(6,*) ' Error in Alloc_obs '  
                end select
                Nt = Ltrot+1-2*Thtrot
                Call Obser_Latt_make(Obs_tau(I),Ns,Nt,No,Filename)
             enddo
          endif
        end Subroutine Alloc_obs

        
!========================================================================
        ! Functions for Global moves.  These move are not implemented in this example.
        Subroutine Global_move(T0_Proposal_ratio,nsigma_old,size_clust)
          
          !>  The input is the field nsigma declared in this module. This routine generates a 
          !>  global update with  and returns the propability  
          !>  T0_Proposal_ratio  =  T0( sigma_out-> sigma_in ) /  T0( sigma_in -> sigma_out)  
          !>   
          Implicit none
          Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio, size_clust
          Integer, dimension(:,:),  allocatable, intent(in)  :: nsigma_old

          T0_Proposal_ratio  = 0.d0

        End Subroutine Global_move
!========================================================================
        Real (Kind=kind(0.d0)) Function Delta_S0_global(Nsigma_old)

          !>  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none 
          
          !> Arguments
          Integer, dimension(:,:), allocatable, intent(IN) :: Nsigma_old

          Delta_S0_global = 0.d0

        end Function Delta_S0_global
        
!========================================================================
        
        Subroutine Obser(GR,Phase,Ntau)
          
          Implicit none
          
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          
          !Local 
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS, ZZ, ZXY
          Integer :: I,J, imj, nf, dec, I1, J1, no_I, no_J,n
          
          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          
          
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Do J = 1,Ndim
                   GRC(I, J, nf) = -GR(J, I, nf)
                Enddo
                GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
             Enddo
          Enddo
          ! GRC(i,j,nf) = < c^{dagger}_{j,nf } c_{j,nf } >

          ! Compute scalar observables. 
          Do I = 1,Size(Obs_scal,1)
             Obs_scal(I)%N         =  Obs_scal(I)%N + 1
             Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo
             

          Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Latt%N
                Do n = 1,N_coord
                   Select Case (Lattice_type)
                   Case ("Square")
                      I1 = I
                      If (n == 1)  J1 = Latt%nnlist(I,1,0)
                      If (n == 2)  J1 = Latt%nnlist(I,0,1)
                   Case ("Honeycomb")
                      I1 = invlist(I,1) 
                      If (n == 1)  J1 = invlist(I,2)
                      If (n == 2)  J1 = invlist(Latt%nnlist(I,1,-1),2)
                      If (n == 3)  J1 = invlist(Latt%nnlist(I,0,-1),2)
                   Case Default
                      stop
                   end Select
                   Zkin = Zkin +  Grc( I1,J1, nf ) + Grc(J1,I1,nf)
                Enddo
             Enddo
          Enddo
          Zkin = -Ham_T*Zkin * dble(N_SUN)
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *ZP* ZS


          dec = 1
          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising" ) then
            dec = 1
          elseif (Model == "Hubbard_MZ") then
            dec = 2
          endif
          Do I = 1,Ndim
             ZPot = ZPot + Grc(i,i,1) * Grc(i,i, dec)
          Enddo
          Zpot = Zpot*ham_U
          Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + Zpot * ZP*ZS


          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Zrho = Zrho + Grc(i,i,nf) 
             enddo
          enddo
          Zrho = Zrho* dble(N_SUN)
          Obs_scal(3)%Obs_vec(1)  =    Obs_scal(3)%Obs_vec(1) + Zrho * ZP*ZS

          Obs_scal(4)%Obs_vec(1)  =    Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*ZP*ZS


          
          ! Compute spin-spin, Green, and den-den correlation functions  !  This is general N_SUN, and  N_FL = 1
          DO I = 1,Size(Obs_eq,1)
             Obs_eq(I)%N        = Obs_eq(I)%N + 1
             Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
          ENDDO
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising" .or. Model == "t_V"  ) then 
             Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   ! Green
                   Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + &
                        &               Z * GRC(I1,J1,1) *  ZP*ZS 
                   ! SpinZ
                   Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + &
                        &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                   ! SpinXY
                   Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + &
                        &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                   ! SpinT
                   Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) + &
                        &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                   ! Den
                   Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J)  +  &
                        &     (    GRC(I1,I1,1) * GRC(J1,J1,1) *Z     + &
                        &          GRC(I1,J1,1) * GR(I1,J1,1 )           &
                        &                                   ) * Z* ZP*ZS
                ENDDO
                Obs_eq(4)%Obs_Latt0(no_I) =  Obs_eq(4)%Obs_Latt0(no_I) +  Z * GRC(I1,I1,1) * ZP * ZS
             ENDDO
          elseif (Model == "Hubbard_Mz" ) Then
             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   ! Green
                   Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + &
                        &                          ( GRC(I1,J1,1) + GRC(I1,J1,2)) *  ZP*ZS 
                   ! SpinZ
                   ZZ =       GRC(I1,J1,1) * GR(I1,J1,1) +  GRC(I1,J1,2) * GR(I1,J1,2)    + &
                        &    (GRC(I1,I1,2) - GRC(I1,I1,1))*(GRC(J1,J1,2) - GRC(J1,J1,1))  
                   Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + ZZ* ZP*ZS
                   ! SpinXY
                   ! c^d_(i,u) c_(i,d) c^d_(j,d) c_(j,u)  +  c^d_(i,d) c_(i,u) c^d_(j,u) c_(j,d)
                   ZXY =  GRC(I1,J1,1) * GR(I1,J1,2) +  GRC(I1,J1,2) * GR(I1,J1,1) 
                   Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + ZXY* ZP*ZS
                   ! SpinT
                   Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) + (2.d0*ZXY + ZZ)*ZP*ZS/3.d0
                   !Den
                   Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) + &
                        & (   GRC(I1,J1,1) * GR(I1,J1,1) +  GRC(I1,J1,2) * GR(I1,J1,2)    + &
                        &   (GRC(I1,I1,2) + GRC(I1,I1,1))*(GRC(J1,J1,2) + GRC(J1,J1,1))     ) * ZP*ZS
                enddo
                Obs_eq(4)%Obs_Latt0(no_I) =  Obs_eq(4)%Obs_Latt0(no_I) +  (GRC(I1,I1,1) + GRC(I1,I1,2)) * ZP*ZS
             enddo
          Endif
                

        end Subroutine Obser
!=====================================================
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE)
          Implicit none
          
          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          
          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS
          Integer :: IMJ, I, J, I1, J1, no_I, no_J

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          If (NT == 0 ) then 
             DO I = 1,Size(Obs_tau,1)
                Obs_tau(I)%N = Obs_tau(I)%N + 1
                Obs_tau(I)%Ave_sign = Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
             ENDDO
          endif
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising" .or. Model == "t_V"  ) then 
             Z =  cmplx(dble(N_SUN),0.d0, kind(0.D0))
             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   ! Green
                   Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        & +  Z * GT0(I1,J1,1) * ZP* ZS
                   
                   ! SpinZ
                   Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        &      - Z*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
                   
                   ! SpinXY
                   Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        &      - Z*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
                   
                   ! Den
                   Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        & + ( Z*Z*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1))*       &
                        &         (cmplx(1.d0,0.d0,kind(0.d0)) - G00(J1,J1,1))  -     &
                        &     Z * GT0(I1,J1,1)*G0T(J1,I1,1)                                ) * ZP * ZS
                Enddo
                Obs_tau(4)%Obs_Latt0(no_I) = Obs_tau(4)%Obs_Latt0(no_I) + &
                     &         Z*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1)) * ZP * ZS
             Enddo
          Elseif ( Model == "Hubbard_Mz"  ) then 
             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   !Green
                   Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        &   +   ( GT0(I1,J1,1) + GT0(I1,J1,2) ) * ZP* ZS

                   !SpinZ
                   Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                       & +  ( &
                       &    (GTT(I1,I1,1) -  GTT(I1,I1,2) ) * ( G00(J1,J1,1)  -  G00(J1,J1,2) )   &
                       &  - (G0T(J1,I1,1) * GT0(I1,J1,1)  +  G0T(J1,I1,2) * GT0(I1,J1,2) )    )*ZP*ZS

                   !SpinXY
                   Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        &  - &
                        &   (G0T(J1,I1,1) * GT0(I1,J1,2)  +  G0T(J1,I1,2) * GT0(I1,J1,1))*ZP*ZS
                   !Den
                   Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        & +  (                                        &  
                        &    (cmplx(2.D0,0.d0,kind(0.d0)) - GTT(I1,I1,1) - GTT(I1,I1,2) ) * &
                        &    (cmplx(2.D0,0.d0,kind(0.d0)) - G00(J1,J1,1) - G00(J1,J1,2) )   &
                        & -  ( G0T(J1,I1,1) * GT0(I1,J1,1) + G0T(J1,I1,2) * GT0(I1,J1,2) )  )*ZP*ZS     

                enddo
             
                Obs_tau(4)%Obs_Latt0(no_I) =  Obs_tau(4)%Obs_Latt0(no_I) + &
                     &       (cmplx(2.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1) - GTT(I1,I1,2)) * ZP * ZS
             Enddo
          Endif
          
        end Subroutine OBSERT

!==========================================================        
        Subroutine  Pr_obs(LTAU)

          Implicit none

          Integer,  Intent(In) ::  Ltau
          
          !Local 
          Integer :: I


          Do I = 1,Size(Obs_scal,1)
             Call  Print_bin_Vec(Obs_scal(I),Group_Comm)
          enddo
          Do I = 1,Size(Obs_eq,1)
             Call  Print_bin_Latt(Obs_eq(I),Latt,dtau,Group_Comm)
          enddo
          If (Ltau  == 1 ) then
             Do I = 1,Size(Obs_tau,1)
                Call  Print_bin_Latt(Obs_tau(I),Latt,dtau,Group_Comm)
             enddo
          endif

        end Subroutine Pr_obs
!===================================================================================           
        Subroutine  Init_obs(Ltau) 

          Implicit none
          Integer, Intent(In) :: Ltau
          
          ! Local 
          Integer :: I

          Do I = 1,Size(Obs_scal,1)
             Call Obser_vec_Init(Obs_scal(I))
          Enddo

          Do I = 1,Size(Obs_eq,1)
             Call Obser_Latt_Init(Obs_eq(I))
          Enddo

          If (Ltau == 1) then
             Do I = 1,Size(Obs_tau,1)
                Call Obser_Latt_Init(Obs_tau(I))
             Enddo
          Endif

        end Subroutine Init_obs
!---------------------------------------------------------------------
        Subroutine  Hamiltonian_set_random_nsigma
          
          ! The user can set the initial configuration
          
          Implicit none
          
          Integer :: I, nt
          
          Do nt = 1,Ltrot
             Do I = 1,Size(OP_V,1)
                nsigma(I,nt)  = 1
                if ( ranf_wrap()  > 0.5D0 ) nsigma(I,nt)  = -1
             enddo
          enddo
          
        end Subroutine Hamiltonian_set_random_nsigma

!---------------------------------------------------------------------
        Subroutine Global_move_tau(T0_Proposal_ratio, S0_ratio, &
             &                     Flip_list, Flip_length,Flip_value,ntau)

!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief 
!> On input: 
!> GR(tau,m) as defined in  Global_tau_mod_PlaceGR and the direction of updating scheme
!> direction=u --> You are visiting the time slices from tau = 1  to tau =Ltrot
!> direction=d --> You are visiting the time slices from tau = Ltrot to tau = 1
!> 
!> On input the field configuration is in the array nsigma.
!> On output: 
!> Flip_list   ::  A list of spins that are to be fliped. Refers to the entires  in OP_V
!> Flip_values ::  The values of the fliped spins
!> Flip_length ::  The number of flips. The first Flip_length entries of Flip_list and Flip_values are relevant
!> S0_ratio          = e^( S_0(sigma_new) ) / e^( S_0(sigma) )
!> T0_Proposal_ratio = T0( sigma_new -> sigma ) /  T0( sigma -> sigma_new)  
!> T0_proposal       = T0 ( sigma -> sigma_new )
!--------------------------------------------------------------------
          
          Implicit none 
          Real (Kind= kind(0.d0)), INTENT(INOUT) :: T0_Proposal_ratio,  S0_ratio
          Integer,    allocatable, INTENT(INOUT) :: Flip_list(:), Flip_value(:)
          Integer, INTENT(INOUT) :: Flip_length
          Integer, INTENT(IN)    :: ntau


          ! Local
          Integer :: n_op, n, ns
          Real (Kind=Kind(0.d0)) :: T0_proposal
          Flip_length = nranf(1)
          
          do n = 1,flip_length 
             n_op = nranf(size(OP_V,1))
             Flip_list(n)  = n_op
             ns            = nsigma(n_op,ntau)
             If ( OP_V(n_op,1)%type == 1 ) then 
                ns = nsigma(n_op,ntau)
                T0_Proposal       =  1.d0 - 1.d0/(1.d0+S0(n_op,ntau)) ! No move prob
                If ( T0_Proposal > Ranf_wrap() ) then
                   T0_Proposal_ratio =  1.d0 / S0(n_op,ntau)
                else
                   T0_Proposal_ratio = 0.d0
                endif
                S0_ratio          =  S0(n_op,ntau)
                Flip_value(n)     = - ns
             else
                Flip_value(n)     = NFLIPL(nsigma(n_op,ntau),nranf(3))
                T0_Proposal_ratio = 1.d0
                S0_ratio          = 1.d0
             endif
          Enddo


        end Subroutine Global_move_tau

      end Module Hamiltonian
