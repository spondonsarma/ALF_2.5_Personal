!  Copyright (C) 2022 The ALF project
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
!> This File is a template for defining new models. 
!> One can define a new model class by copying this file, replacing alle occurences
!> of ##NAME## by the Hamiltonian name, populating the subroutines below as needed
!> adding the Hamiltonian name to the file Prog/Hamiltonians.list.

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
!>@param [public] WF_L
!> \verbatim Type (WaveFunction), dimension(:),   allocatable
!> Left trial wave function.  \endverbatim
!>
!> @param [public] WF_R
!> \verbatim Type (WaveFunction), dimension(:),   allocatable
!> Right trial wave function.   For both wave functions the index runs over the flavor index. \endverbatim
!>
!> @param [public]  nsigma
!> \verbatim Type(Fields)
!> Contains all auxiliary fields in the variable f(:,:). The first index runs through the operator
!> sequence. The second through the time slices.   \endverbatim
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
!> symmetrically. If Symm is true, the propagation reads:
!> \f$ \prod_{\tau} \; \;  \prod_{n=N_T}^{1}e^{T_n/2} \prod_{n=1}^{N_V}e^{V_n(\tau)}  \prod_{n=1}^{N_T}e^{T_n/2}  \f$
!>
!>
!> You still have to add some docu for the other private variables in this module.
!>
!--------------------------------------------------------------------

    submodule (Hamiltonian_main) ham_Kondo_dim_mismatch_smod

      Use Operator_mod
      Use WaveFunction_mod
      Use Lattices_v3
      Use MyMats
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      Use Observables
      Use Fields_mod
      Use Predefined_Hoppings
      Use Predefined_Int
      Use LRC_Mod

      Implicit none
      
      type, extends(ham_base) :: ham_Kondo_dim_mismatch
      contains
        ! Set Hamiltonian-specific procedures
        procedure, nopass :: Ham_Set
        procedure, nopass :: Alloc_obs
        procedure, nopass :: Obser
        procedure, nopass :: ObserT
        procedure, nopass :: Ham_Latt
        procedure, nopass :: Ham_Hop
#ifdef HDF5
        procedure, nopass :: write_parameters_hdf5
#endif
      end type ham_Kondo_dim_mismatch

      !#PARAMETERS START# VAR_lattice
      Character (len=64) :: Model = ''  ! Value not relevant
      Character (len=64) :: Lattice_type = 'Square'  ! Possible values: 'Bilayer_square', 'Bilayer_honeycomb'
      Integer            :: L1 = 6   ! Length in direction a_1
      Integer            :: L2 = 6   ! Length in direction a_2
      !#PARAMETERS END#

      !#PARAMETERS START# VAR_Model_Generic
      !Integer              :: N_SUN        = 2        ! Number of colors
      !Integer              :: N_FL         = 1        ! Number of flavors
      real(Kind=Kind(0.d0)) :: Phi_X        = 0.d0     ! Twist along the L_1 direction, in units of the flux quanta
      real(Kind=Kind(0.d0)) :: Phi_Y        = 0.d0     ! Twist along the L_2 direction, in units of the flux quanta
      logical               :: Bulk         = .true.   ! Twist as a vector potential (.T.), or at the boundary (.F.)
      Integer               :: N_Phi        = 0        ! Total number of flux quanta traversing the lattice
      real(Kind=Kind(0.d0)) :: Dtau         = 0.1d0    ! Thereby Ltrot=Beta/dtau
      real(Kind=Kind(0.d0)) :: Beta         = 5.d0     ! Inverse temperature
      logical               :: Checkerboard = .true.   ! Whether checkerboard decomposition is used
      !logical              :: Symm         = .true.   ! Whether symmetrization takes place
      !logical              :: Projector    = .false.  ! Whether the projective algorithm is used
      real(Kind=Kind(0.d0)) :: Theta        = 10.d0    ! Projection parameter
      !#PARAMETERS END#

      !#PARAMETERS START# VAR_Kondo_dim_mismatch
      real(Kind=Kind(0.d0)) :: ham_T    = 1.d0  ! Hopping parameter
      real(Kind=Kind(0.d0)) :: Ham_chem = 0.d0  ! Chemical potential
      real(Kind=Kind(0.d0)) :: Ham_JK   = 2.d0  ! Kondo Coupling  J
      real(Kind=Kind(0.d0)) :: Ham_Jh   = 1.d0  ! Heisenberg  coupling
      real(Kind=Kind(0.d0)) :: Ham_U    = 1.d0  ! Hubbard
      !#PARAMETERS END#

      Type (Lattice),       target :: Latt_c
      Type (Unit_cell),     target :: Latt_unit_c
      Type (Hopping_Matrix_type), Allocatable :: Hopping_Matrix(:)
      Integer, allocatable :: List_c(:,:), Invlist_c(:,:)  

      Type (Lattice),       target :: Latt_f
      Type (Unit_cell),     target :: Latt_unit_f
      Integer, allocatable :: List_f(:,:), Invlist_f(:,:)  

    contains
      
      module Subroutine Ham_Alloc_Kondo_dim_mismatch
        allocate(ham_Kondo_dim_mismatch::ham)
      end Subroutine Ham_Alloc_Kondo_dim_mismatch

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_Kondo_dim_mismatch_read_write_parameters.F90"

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the Hamiltonian
!--------------------------------------------------------------------
      Subroutine Ham_Set

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif
          Implicit none

          integer                :: ierr, nf, unit_info
          Character (len=64)     :: file_info


#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif

!         !From dynamically generated file "Hamiltonian_Kondo_dim_mismatch_read_write_parameters.F90"
          call read_parameters()
 
          Ltrot = nint(beta/dtau)
          if (Projector) Thtrot = nint(theta/dtau)
          Ltrot = Ltrot+2*Thtrot
 
!         ! Setup the Bravais lattice
          call Ham_Latt
! 
!           ! Setup the hopping / single-particle part
          call Ham_Hop
! 
!           ! Setup the interaction.
          call Ham_V
! 
!           ! Setup the trival wave function, in case of a projector approach
!           if (Projector) Call Ham_Trial()

#ifdef MPI
          If (Irank_g == 0) then
#endif
             File_info = "info"
#if defined(TEMPERING)
             write(File_info,'(A,I0,A)') "Temp_",igroup,"/info"
#endif
             Open(newunit=unit_info, file=file_info, status="unknown", position="append")
             Write(unit_info,*) '====================================='
             Write(unit_info,*) 'Model is      :  Kondo_dim_mismatch'
             Write(unit_info,*) 'Lattice is    : ', Lattice_type
             Write(unit_info,*) '# c-orbs      : ', Latt_c%N 
             Write(unit_info,*) '# f-orbs      : ', Latt_f%N
             Write(unit_info,*) '# orbs        : ', Ndim
             Write(unit_info,*) '# t           : ', Ham_t
             Write(unit_info,*) '# J_k         : ', Ham_Jk
             Write(unit_info,*) '# J_h         : ', Ham_Jh
             Write(unit_info,*) '# U           : ', Ham_U
             Write(unit_info,*) '# Flux        : ', N_Phi
             Write(unit_info,*) 'Checkerboard  : ', Checkerboard
             Write(unit_info,*) 'Symm. decomp  : ', Symm
             if (Projector) then
                Write(unit_info,*) 'Projective version'
                Write(unit_info,*) 'Theta         : ', Theta
                Write(unit_info,*) 'Tau_max       : ', beta
             else
                Write(unit_info,*) 'Finite temperture version'
                Write(unit_info,*) 'Beta          : ', Beta
             endif
             Write(unit_info,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
             Write(unit_info,*) 'N_SUN         : ',   N_SUN
             Write(unit_info,*) 'N_FL          : ', N_FL
             if (Projector) then
                Do nf = 1,N_FL
                   Write(unit_info,*) 'Degen of right trial wave function: ', WF_R(nf)%Degen
                   Write(unit_info,*) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
                enddo
             endif
             Close(unit_info)
#ifdef MPI
          Endif
#endif
        end Subroutine Ham_Set

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Specifies the lattice
!> @details
!--------------------------------------------------------------------
        Subroutine  Ham_Latt

          Implicit  none
          
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2), ic_p(2)
          Integer ::  nc, I, no, Ic

          Latt_Unit_c%Norb      = 1
          Allocate (Latt_unit_c%Orb_pos_p(1,2))
          Latt_Unit_c%Orb_pos_p(1,:) = 0.d0
          a1_p(1) =  1.0  ; a1_p(2) =  0.d0
          a2_p(1) =  0.0  ; a2_p(2) =  1.d0
          L1_p    =  dble(L1)*a1_p
          L2_p    =  dble(L2)*a2_p
          Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt_c )
          
          

          Latt_Unit_f%Norb      = 2
          Allocate (Latt_unit_f%Orb_pos_p(1,2))
          Latt_Unit_c%Orb_pos_p(1,:) = 0.d0
          a1_p(1) =  1.0  ; a1_p(2) =  0.d0
          a2_p(1) =  0.0  ; a2_p(2) =  1.d0
          L1_p    =  dble(L1)*a1_p
          L2_p    =  a2_p
          Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt_f )
          
          Ndim  =  Latt_c%N  + Latt_f%N

          !  Setup  the lists
          Allocate (List_c(Latt_c%N*Latt_Unit_c%Norb,2), Invlist_c(Latt_c%N,Latt_Unit_c%Norb))
          Allocate (Invlist_f(Latt_f%N,Latt_Unit_f%Norb))
          nc = 0
          Do I = 1,Latt_c%N
             Do no = 1,Latt_Unit_c%Norb
                nc = nc + 1
                List_c   (nc,1) = I
                List_c   (nc,2) = no
                Invlist_c(I,no) = nc
             Enddo
          Enddo
          Do I = 1,Latt_f%N
             nc   = nc + 1
             ic_p = dble(Latt_f%list(I,1))*Latt_f%a1_p + dble(Latt_f%list(I,2))*Latt_f%a2_p
             Ic   = Inv_R(ic_p,Latt_c)
             Invlist_f(I,1) = nc                  !  f-orbital
             Invlist_f(I,2) = Invlist_c(Ic,1)     !  c-orbital
          enddo

          !Testing
          Do I =  1,Latt_f%N
             ic_p = dble(Latt_f%list(I,1))*Latt_f%a1_p + dble(Latt_f%list(I,2))*Latt_f%a2_p
             Write(6,"(I4,2x,F14.7,2x,F14.7)") I, ic_p(1), ic_p(2)
             Ic   = Inv_R(ic_p,Latt_c)
             ic_p = dble(Latt_c%list(Ic,1))*Latt_c%a1_p + dble(Latt_c%list(Ic,2))*Latt_c%a2_p
             Write(6,"(I4,2x,F14.7,2x,F14.7)") Ic, ic_p(1), ic_p(2)
             Write(6,*)
          enddo
          
        end Subroutine Ham_Latt
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the  hopping
!> @details
!--------------------------------------------------------------------
        Subroutine  Ham_Hop

          Implicit  none

          Real (Kind=Kind(0.d0) ), allocatable :: Ham_T_vec(:),  Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:)
          Integer, allocatable ::   N_Phi_vec(:)


          Select case  (Lattice_type)
          Case("Square")
             N_phi  =  0
          Case ("Pi_flux")
             N_phi  =  Latt_c%N/2
          case default
             Write(error_unit,*) 'Lattice  has  to be  set to Pi_flux or  Square'
             CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
          end Select


          Allocate (Ham_T_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL), N_Phi_vec(N_FL) )

          ! Here we consider no N_FL  dependence of the hopping parameters.
          Ham_T_vec      = Ham_T
          Ham_Chem_vec   = Ham_Chem
          Phi_X_vec      = Phi_X
          Phi_Y_vec      = Phi_Y
          N_Phi_vec      = N_Phi

          Call  Set_Default_hopping_parameters_square(Hopping_Matrix,Ham_T_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
               &                                      Bulk, N_Phi_vec, N_FL, List_c, Invlist_c, Latt_c, Latt_unit_c )
          
          Call  Predefined_Hoppings_set_OPT(Hopping_Matrix,List_c,Invlist_c,Latt_c,  Latt_unit_c,  Dtau, Checkerboard, Symm, OP_T )

          Deallocate (Ham_T_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, N_Phi_vec )

        end Subroutine Ham_Hop
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the interaction.
!> @details
!> The  Hubbard  interaction  commutes  with the   Heisenberg and
!> Kondo  interactions, such  that we  can pull it out.
!> For  the  other  terms  we  use  a  symmettic    decompostion. 
!> e^{-\Delta \tau H_J^{even}/2} e^{-\Delta \tau H_J^{odd}/2}  e^{-Delta \tau  H_JK}
!> e^{-\Delta \tau H_J^{odd}/2} e^{-\Delta \tau H_J^{even}/2} e^{-\Delta \tau H_U}
!--------------------------------------------------------------------
        Subroutine  Ham_V

          Implicit  none

          Integer :: N_op_U ,  N_op_K, N_op_H,  N_op,  n,  nf, nc
          Integer :: I, J 
          
          N_op_U   =    Latt_f%N  !  Hubbard
          N_op_K   =    Latt_f%N  !  Kondo
          N_op_H   =  2*Latt_f%N  !  Heisenberg

          N_op =  N_op_U +  N_op_K +   N_op_H
          Allocate(Op_V(N_op,N_FL))
          Do  nf =  1,N_FL
             nc = 0
             ! Hubbard
             Do  n = 1,Latt_f%N
                nc = nc + 1
                I = Invlist_f(n,1)  !  f-orbital 
                Call Predefined_Int_U_SUN( OP_V(nc,nf), I, N_SUN, DTAU, Ham_U  )
             enddo
             ! Heisenberg
             Do n  =  1,Latt_f%N
                nc = nc + 1
                I = Invlist_f( n                   ,1)
                J = Invlist_f( latt_f%nnlist(n,1,0),1)
                Call Predefined_Int_V_SUN( OP_V(nc,nf), I, J, N_SUN, DTAU, Ham_Jh/4.d0  ) 
             enddo
             ! Kondo
             Do n  =  1,Latt_f%N
                nc = nc + 1
                I = Invlist_f( n                   ,1)
                J = Invlist_f( n                   ,2)
                Call Predefined_Int_V_SUN( OP_V(nc,nf), I, J, N_SUN, DTAU, Ham_Jk/2.d0  ) 
             enddo
             ! Heisenberg
             Do n  =  Latt_f%N,1, -1
                nc = nc + 1
                I = Invlist_f( n                   ,1)
                J = Invlist_f( latt_f%nnlist(n,1,0),1)
                Call Predefined_Int_V_SUN( OP_V(nc,nf), I, J, N_SUN, DTAU, Ham_Jh/4.d0  ) 
             enddo
          enddo

          
        End Subroutine Ham_V
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Specifiy the equal time and time displaced observables
!> @details
!--------------------------------------------------------------------
        Subroutine  Alloc_obs(Ltau)

          Implicit none
          !>  Ltau=1 if time displaced correlations are considered.
          Integer, Intent(In) :: Ltau
          Integer    ::  i, N, Nt
          Character (len=64) ::  Filename
          Character (len=2)  ::  Channel


!         Scalar observables
          Allocate ( Obs_scal(1) )
          Do I = 1,Size(Obs_scal,1)
             select case (I)
             case (1)
                N = 1;   Filename = "Pot"
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
          enddo
          
!         ! Equal time correlators
          Allocate ( Obs_eq(1) )
          Do I = 1,Size(Obs_eq,1)
             select case (I)
             case (1)
                Filename = "Spin"
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
             Nt = 1
             Channel = '--'
             Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt_f, Latt_unit_c, Channel, dtau)
          enddo
! 
!           If (Ltau == 1) then
!             ! Time-displaced correlators
!             Allocate ( Obs_tau(3) )
!             Do I = 1,Size(Obs_tau,1)
!               select case (I)
!               case (1)
!                 Channel = 'P' ; Filename = "Green"
!               case (2)
!                 Channel = 'PH'; Filename = "SpinZ"
!               case (3)
!                 Channel = 'PH'; Filename = "Den"
!               case default
!                 Write(6,*) ' Error in Alloc_obs '
!               end select
!               Nt = Ltrot+1-2*Thtrot
!               If(Projector) Channel = 'T0'
!               Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
!             enddo
!           endif

        End Subroutine Alloc_obs

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Computes equal time observables
!> @details
!> @param [IN] Gr   Complex(:,:,:)
!> \verbatim
!>  Green function: Gr(I,J,nf) = <c_{I,nf } c^{dagger}_{J,nf } > on time slice ntau
!> \endverbatim
!> @param [IN] Phase   Complex
!> \verbatim
!>  Phase
!> \endverbatim
!> @param [IN] Ntau Integer
!> \verbatim
!>  Time slice
!> \endverbatim
!-------------------------------------------------------------------
        subroutine Obser(GR,Phase,Ntau, Mc_step_weight)

          Use Predefined_Obs

          Implicit none

          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer,                   INTENT(IN) :: Ntau
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

          !Local
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)) :: ZP, ZS
          Complex (Kind=Kind(0.d0)) :: Z_pot, Z
          Integer :: I, I1,  J, J1,  nf, n, imj
          ! Add local variables as needed

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))

          ZS = ZS*Mc_step_weight
          
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Do J = 1,Ndim
                   GRC(I, J, nf) = -GR(J, I, nf)
                Enddo
                GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
             Enddo
          Enddo
          
          ! GRC(i,j,nf) = < c^{dagger}_{i,nf } c_{j,nf } >
          ! Compute scalar observables.
          ! Compute scalar observables.
          Do I = 1,Size(Obs_scal,1)
             Obs_scal(I)%N         =  Obs_scal(I)%N + 1
             Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo

          Z_pot  = cmplx(0.d0,0.d0, kind(0.D0))
          do I  =  1, Latt_f%N
             I1  =  invlist_f(I,1)  !  f-orbital
             Z_pot  = Z_pot  +  Grc(I1,I1,1)*Grc(I1,I1,1)
          enddo
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Z_pot *ZP* ZS
          
          ! Compute equal-time correlations
          Do I = 1,Size(Obs_eq,1)
             Obs_eq(I)%N         =  Obs_eq(I)%N + 1
             Obs_eq(I)%Ave_sign  =  Obs_eq(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo
          
          Do I  = 1, Latt_f%N
             I1  =   Invlist_f(I,1) !  f-orbital 
             Do J = 1,  Latt_f%N
                J1  =   Invlist_f(J,1) !  f-orbital 
                imj = latt_f%imj(I,J)
                Z =  GRC(I1,J1,1) * GR(I1,J1,1) * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
                Obs_eq(1)%Obs_Latt(imj,1,1,1) =  Obs_eq(1)%Obs_Latt(imj,1,1,1) + Z*ZP*ZS
             Enddo
          Enddo
          
          
        end Subroutine Obser


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Computes time displaced  observables
!> @details
!> @param [IN] NT, Integer
!> \verbatim
!>  Imaginary time
!> \endverbatim
!> @param [IN] GT0, GTT, G00, GTT,  Complex(:,:,:)
!> \verbatim
!>  Green functions:
!>  GT0(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(0  )>
!>  G0T(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(tau)>
!>  G00(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(0  )>
!>  GTT(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(tau)>
!> \endverbatim
!> @param [IN] Phase   Complex
!> \verbatim
!>  Phase
!> \endverbatim
!-------------------------------------------------------------------
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE,  Mc_step_weight)

          Use Predefined_Obs

          Implicit none

          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight
          
          !Locals
          Complex (Kind=Kind(0.d0)) :: ZP, ZS
          ! Add local variables as needed

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          ZS = ZS * Mc_step_weight

          ! Compute observables

        end Subroutine OBSERT
        
      end submodule ham_Kondo_dim_mismatch_smod
