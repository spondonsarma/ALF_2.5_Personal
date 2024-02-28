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

    submodule (Hamiltonian_main) ham_Kondo_Impurity_Wilson_smod

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
      Use Predefined_Obs
      Use LRC_Mod
      use runtime_error_mod

      Implicit none
      
      type, extends(ham_base) :: ham_Kondo_Impurity_Wilson
      contains
        ! Set Hamiltonian-specific procedures
        procedure, nopass :: Ham_Set
        procedure, nopass :: Alloc_obs
        procedure, nopass :: Obser
        procedure, nopass :: ObserT
        procedure, nopass :: weight_reconstruction
        procedure, nopass :: GR_reconstruction
        procedure, nopass :: GRT_reconstruction
        ! procedure, nopass :: ##PROCEDURE_NAME##  ! Some other procedure defined in ham_base
#ifdef HDF5
        procedure, nopass :: write_parameters_hdf5
#endif
      end type ham_Kondo_Impurity_Wilson

      Type (Lattice),       target :: Latt
      Type (Unit_cell),     target :: Latt_unit
      Type (Hopping_Matrix_type), Allocatable :: Hopping_Matrix(:)
      Integer, allocatable :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell

      
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

      !#PARAMETERS START# VAR_Kondo_Impurity_Wilson
      Integer                  :: L_Bath = 12        !   Discretization
      Real(Kind=Kind(0.d0))    :: Ham_W = 1.d0      !   Bandwidth
      Real(Kind=Kind(0.d0))    :: Lambda = 1.d0 !   For  log  discretization
      Integer                  :: Two_S  = 1    !   2spin of  impurity
      Real(Kind=Kind(0.d0))    :: Ham_JK = 1.d0     !   Kondo
      Real(Kind=Kind(0.d0))    :: Ham_D  = 1.d0     !   Easyplane
      Real(Kind=Kind(0.d0))    :: Ham_Jh = -1.d0    !   Heisenberg
      Real(Kind=Kind(0.d0))    :: Ham_U  =  1.d0    !    Hubbard
      !#PARAMETERS END#

      real (Kind=Kind(0.d0)),    allocatable ::  eps(:), delta_eps(:), g(:)
      INTEGER :: nf_calc, nf_reconst
      Integer, allocatable   ::  Prefactor(:)
      Logical  :: Particle_hole

    contains

      module Subroutine Ham_Alloc_Kondo_Impurity_Wilson
        allocate(ham_Kondo_Impurity_Wilson::ham)
      end Subroutine Ham_Alloc_Kondo_Impurity_Wilson

      ! Dynamically generated on compile time from parameters list.
      ! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_Kondo_Impurity_Wilson_read_write_parameters.F90"
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

          integer                :: ierr, nf, unit_info, i
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

!       ! From dynamically generated file "Hamiltonian_Kondo_Impurity_Wilson_read_write_parameters.F90"
        call read_parameters()
        Ltrot = nint(beta/dtau)
        Thtrot = 0
        if (Projector) Thtrot = nint(theta/dtau)
        Ltrot = Ltrot+2*Thtrot
        NDim   =  L_Bath + Two_S
        ! Setup  lists  for log   discretization (This plays  the  role of the lattice)
        ! Bath sites   1 .. L_Bath,   Impurity  sites   L_Bath + 1,  L + 2S  
        call Ham_Bath
        ! Setup diagonal hopping =  site  dependent  energy
        call Ham_Hop
        ! Setup the interaction.
        call Ham_V
        ! Allocate a lattice of one  unit  cell  with  
        ! Norb  =  Two_S
        Call  Ham_Latt
        if  (N_FL == 2)  then
            !Setup  prefactor  for falvor  symmetry
           Call  Set_Prefactor(Particle_hole)
         !  Particle_hole =.false.
           If  (Particle_hole)  then
              allocate(Calc_Fl(N_FL))
              nf_calc=2
              nf_reconst=1
              Calc_Fl(nf_calc)=.True.
              Calc_Fl(nf_reconst)=.False.
           endif
        endif
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
          Write(unit_info,*) 'Model is      : Kondo_Impurity_Wilson '
          Write(unit_info,*) '2S            : ',Two_S
          Write(unit_info,*) 'L Bath        : ',L_Bath
          Write(unit_info,*) 'Ndim          : ',Ndim
          Write(unit_info,*) 'W             : ',Ham_W
          Write(unit_info,*) 'Jk            : ',Ham_JK
          Write(unit_info,*) 'U             : ',Ham_U
          Write(unit_info,*) 'Lambda        : ',Lambda
          Write(unit_info,*) 'Beta          : ', Beta
          Write(unit_info,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
          Write(unit_info,*) 'N_SUN         : ', N_SUN
          Write(unit_info,*) 'N_FL          : ', N_FL
!         if (Projector) then
!            Do nf = 1,N_FL
!               Write(unit_info,*) 'Degen of right trial wave function: ', WF_R(nf)%Degen
!               Write(unit_info,*) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
!            enddo
!         endif
          Close(unit_info)
#ifdef MPI
        Endif
#endif

#ifdef MPI
        If (Irank_g == 0 ) then
#endif
          Open (Unit=10,file="Discretizazion",status="unknown")
          Write(10,*) 'Disretization'
          Do i = 1,L_Bath
            Write(10,*)  i, eps(i), delta_eps(i)
          enddo
          close(10)
#ifdef MPI
        endif
#endif
       

      end Subroutine Ham_Set
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>      Sets the   required  tables  for   the  discretization of the  density of  states
!> @details
!--------------------------------------------------------------------
      Subroutine  Ham_Bath
        Implicit none

        Integer ::  I
        Real (Kind=Kind(0.d0)) :: Delta
        Logical :: Log_scale  = .true. , Lin_scale = .false.

        !Here you can allocate the energies
        allocate ( eps(L_Bath), delta_eps(L_Bath), g(L_Bath) )
        Do I = 1,L_Bath
           g(i)  = real(L_Bath,kind(0.d0)) / ( Ham_W * 2.D0 )
        Enddo
        if (Log_scale) then
          Do I = 1,L_Bath/2
            eps(I) =  - (lambda**( 1 -i) ) * Ham_W/2.d0
          Enddo
          Do I = L_Bath/2+1,L_Bath
            eps(I) =   (lambda**(i-L_Bath)) * Ham_W/2.d0
          Enddo
        elseif (Lin_scale) then
          Delta = Ham_W/real(L_Bath,kind(0.d0))
          Do I = 1,L_Bath/2
            eps(I) =  -Ham_W/2.d0   + real(i-1,kind(0.d0))*Delta
          Enddo
          Do I = L_Bath/2+1,L_Bath
            eps(I) =   Delta*real( I - L_Bath/2,Kind(0.d0))
          Enddo
        else
          Write(6,*) 'No valid scale'
          Call Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
        endif
        Do I = 1,L_Bath/2-1
           delta_eps(I) = abs(eps(I+1) - eps(I))
        Enddo
        delta_eps(L_Bath/2)   = abs(eps(L_Bath/2))
        delta_eps(L_Bath/2+1) = abs(eps(L_Bath/2))
        Do I = L_Bath/2+2,L_Bath
           delta_eps(I) = abs(eps(I-1) - eps(I))
        Enddo
      end   subroutine Ham_Bath

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the  Lattice. 
!--------------------------------------------------------------------
      Subroutine Ham_Latt

        Implicit none
        Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
        
        
        Latt_unit%Norb    = Two_S
        Latt_unit%N_coord = 1
        allocate(Latt_unit%Orb_pos_p(Latt_unit%Norb,2))
        Latt_unit%Orb_pos_p(1, :) = [0.d0, 0.d0]

        a1_p(1) =  1.0  ; a1_p(2) =  0.d0
        a2_p(1) =  0.0  ; a2_p(2) =  1.d0
        L1_p    =  a1_p
        L2_p    =  a2_p
        Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )


       
        
      end Subroutine Ham_Latt
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>       Specifies  the  hoping.  In our  case  this is just  a  site  dependent  diagonal term. 
!> @details
!--------------------------------------------------------------------
      Subroutine Ham_hop
        Implicit none

        Integer :: I , nf

        allocate(Op_T(L_Bath,N_FL))
        do nf = 1,N_Fl
          Do I = 1,L_Bath
             Call Op_make(Op_T(I,nf),1)
            Op_T(I,nf)%O(1,1) = cmplx(eps(i)*delta_eps(i)*g(i), 0.d0, kind(0.D0))
            Op_T(I,nf)%P(1) =  I
            Op_T(I,nf)%g      = -Dtau
            Op_T(I,nf)%alpha  = cmplx( 0.d0, 0.d0, kind(0.D0) )
            Call Op_set( Op_T(I,nf) )
          enddo
        enddo

      end  subroutine Ham_hop
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!>       Sets  the  interaction
!> @details
!--------------------------------------------------------------------
      Subroutine Ham_V

        Implicit none

        
        Integer ::  N_op,  nc,  Nf, n, I, I1, Spin
        Real  (Kind=Kind(0.d0)) :: Sign_D
        N_op  =  Two_S   +  Two_S  + (Two_S - 1)    +  (Two_S -1 )
        !        J_K        U        J_h                  Easy axis       
        Allocate (Op_V(N_op,N_Fl))
        do nf =  1, N_fl
          Spin = 1
          if  (nf == 2) Spin = -1
          nc = 0 
          do n = 1, Two_S  ! *******  Kondo   *******
            nc = nc + 1
            Call Op_make(Op_V(nc,nf),L_Bath+1)    
            Op_V(nc,nf)%P(L_Bath+1) = L_Bath + n
            Do I = 1, L_Bath
              Op_V(nc,nf)%P(I ) = I
              Op_V(nc,nf)%O(I,L_Bath+1) = cmplx(g(i)*Delta_eps(i)  ,0.d0, kind(0.D0))
              Op_V(nc,nf)%O(L_Bath+1,I) = cmplx(g(i)*Delta_eps(i)  ,0.d0, kind(0.D0))
            Enddo
            Op_V(nc,nf)%alpha  = cmplx(0.d0  ,0.d0, kind(0.d0) )
            Op_V(nc,nf)%g      = cmplx(sqrt(2.d0*Dtau*Ham_Jk/(4.d0*Real(L_Bath,Kind(0.d0)))),0.d0, kind(0.D0))
            Op_V(nc,nf)%type   = 2
            Call Op_set( Op_V(nc,nf) )
          enddo

          Do n =  1, Two_S ! Hubbard
            nc = nc + 1
            Call Op_make(Op_V (nc,nf),1) 
            Op_V(nc,nf)%P(1) =  L_Bath + n
            Op_V(nc,nf)%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
            Op_V(nc,nf)%alpha  = cmplx(-0.5d0,0.d0, kind(0.D0))
            Op_V(nc,nf)%g      = SQRT(CMPLX(-DTAU*ham_U/2.d0, 0.D0, kind(0.D0)))
            Op_V(nc,nf)%type   = 2
            Call Op_set( Op_V(nc,nf) )
          enddo
          do n = 1,Two_S -1   !   Ferromagnetic  Heisenberg.
            nc = nc + 1
            Call Op_make(Op_V(nc,nf),2  )    
            I  = L_Bath + n
            I1 = L_Bath + n + 1
            Op_V(nc,nf)%P(1)   = I
            Op_V(nc,nf)%P(2)   = I1
            Op_V(nc,nf)%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.D0)) 
            Op_V(nc,nf)%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
            Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_Jh/4.d0, 0.D0, kind(0.D0))) 
            Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
            Op_V(nc,nf)%type   = 2
            Call Op_set( Op_V(nc,nf)  )
          enddo
          

          do n = 1,Two_S -1   !  Easy   axis.  
            nc = nc + 1
            Call Op_make(Op_V (nc,nf),2) 
            I  = L_Bath+n
            I1 = L_Bath+n+1    

            Sign_D = 1.d0
            If  (Abs(Ham_D) > 10.D-10)   Sign_D = Ham_D/Abs(Ham_D)

            Op_V(nc,nf)%P(1)   =  I
            Op_V(nc,nf)%P(2)   =  I1
            Op_V(nc,nf)%O(1,1) =  cmplx(dble(Spin)  ,0.d0, kind(0.D0))
            Op_V(nc,nf)%O(2,2) =  cmplx(-dble(Spin)*Sign_D  ,0.d0, kind(0.D0))
            Op_V(nc,nf)%alpha  =  cmplx(0.d0, 0.d0, kind(0.D0))
            Op_V(nc,nf)%g      =  SQRT(CMPLX(DTAU*Abs(Ham_D)/8.d0, 0.D0, kind(0.D0))) 
            Op_V(nc,nf)%type   =  2
            Call Op_set( Op_V(nc,nf)  )
            enddo
        enddo
        Write(6,*)  'Total  # of  local  operators  for  interaction ',  nc, N_op
      end Subroutine  Ham_V

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Specifiy the equal time and time displaced observables
!> @details
!--------------------------------------------------------------------
Subroutine Alloc_obs(Ltau)

  Implicit none
  !>  Ltau=1 if time displaced correlations are considered.
  Integer, Intent(In) :: Ltau
  Integer    ::  i, N, Nt
  Character (len=64) ::  Filename
  Character (len=2)  ::  Channel

  ! Scalar observables
  Allocate ( Obs_scal(1) )
  Do I = 1,Size(Obs_scal,1)
    select case (I)
    case (1)
      N = 1;   Filename = "Pot"
    case default
      Write(6,*) ' Error in Alloc_obs '
      Call Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
    end select
    Call Obser_Vec_make(Obs_scal(I),N,Filename)
  enddo

  ! Scalar correlators
  ! Allocate ( Obs_eq(1) )
  ! Do I = 1,Size(Obs_eq,1)
  !   select case (I)
  !   case (1)
  !     Filename = "SpinTOT"
  !   case default
  !     Write(6,*) ' Error in Alloc_obs '
  !     Call Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
  !   end select
  !   Nt = 1
  !   Channel = '--'
  !   Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
  ! enddo

          
!         Equal time correlators
          If  (N_SUN == 1)  then
             Allocate ( Obs_eq(3) )
             Do I = 1,Size(Obs_eq,1)
                select case (I)
                case (1)
                   Filename = "SpinZ"
                case (2)
                   Filename = "SpinXY"
                case (3)
                   Filename = "SpinT"
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = 1
                Channel = '--'
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             enddo
          else
             Allocate ( Obs_eq(1) )
             Do I = 1,Size(Obs_eq,1)
                select case (I)
                case (1)
                   Filename = "SpinZ"
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = 1
                Channel = '--'
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             enddo
          endif
  If (Ltau == 1) then
    ! Time-displaced correlators
    If  (N_FL == 2)  then
      Allocate ( Obs_tau(4) )
    Else
      Allocate ( Obs_tau(2) )
    endif

    If  (N_FL == 2)  then
      Do I = 1,Size(Obs_tau,1)
        select case (I)
        case (1)
          Channel = 'PH' ; Filename = "SpinZ"
        case (2)
          Channel = 'PH' ; Filename = "SpinXY"
        case (3)
          Channel = 'PH' ; Filename = "SpinT"
        case (4)
          Channel = 'P'; Filename = "Psi"
        case default
          Write(6,*) ' Error in Alloc_obs '
          Call Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
        end select
        Nt = Ltrot+1-2*Thtrot
        If(Projector) Channel = 'T0'
        Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
      enddo
   

  else
    Do I = 1,Size(Obs_tau,1)
      select case (I)
      case (1)
        Channel = 'PH'; Filename = "SpinZ"
      case (2)
        Channel = 'P'; Filename = "Psi"
      case default
        Write(6,*) ' Error in Alloc_obs '
        Call Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
      end select
      Nt = Ltrot+1-2*Thtrot
      If(Projector) Channel = 'T0'
      Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
    enddo
  endif
 endif

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
          Complex (Kind=Kind(0.d0)) :: ZP, ZS,  ZPOT, ZZ, ZXY
          Integer :: I, J, nf,n,m,f
          Complex (Kind=Kind(0.d0)):: Total_Spin!, Total_Spin_Squared
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
          Do I = 1,Size(Obs_scal,1)
             Obs_scal(I)%N         =  Obs_scal(I)%N + 1
             Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo
         f=1


         ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
         Do I = L_Bath+f, L_Bath + Two_S
         If (N_FL==1) then
            ZPot = ZPot + GRC(I,I,1) * GRC(I,I,1)
         else
            ZPot= Zpot + GRC (I,I,1)* GRC (I,I,2)    
         endif 
         enddo
         Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + ZPot *ZP* ZS
    
          
          !   Obs_eq(1)%N         =  Obs_eq(1)%N + 1
         !    Obs_eq(1)%Ave_sign  =  Obs_eq(1)%Ave_sign + Real(ZS,kind(0.d0))
         

       !   Do n = 1,Two_S-1
        !  I=L_Bath+n
         ! J=L_Bath+n+1    
        !  Total_Spin=Total_Spin+ (GRC(I, I, 1)) * (GRC(J, J, 1))
        !  Obs_eq(1)%Obs_Latt(1,1,1,1)  = Obs_eq(1)%Obs_Latt(1,1,1,1) +2.d0*Total_Spin*ZS*ZP
        !  enddo
             
        

          !Total_Spin_Squared =Total_Spin**2

         If (N_SUN == 2 )  then    
          
         Do I = 1, Size(Obs_eq, 1)
            Obs_eq(I)%N         = Obs_eq(I)%N + 1
            Obs_eq(I)%Ave_sign  = Obs_eq(I)%Ave_sign + Real(ZS, kind(0.d0)) 
         Enddo

         Do n = 1, Two_S
            I = L_Bath + n
          Do m = 1, Two_S
            J = L_Bath + m
            ZZ = 2.d0 * GRC(I, J, 1) * GR(I, J, 1)
            Obs_eq(1)%Obs_Latt(1, 1, n, m) = Obs_eq(1)%Obs_Latt(1, 1, n, m) + ZZ * ZP * ZS
            ! Display the values of I and J
            !Write(6,*)  'I,J= ',  I, J
          Enddo
         Enddo

 
         else
          
          Do I = 1,Size(Obs_eq,1)
               Obs_eq(I)%N         =  Obs_eq(I)%N + 1
               Obs_eq(I)%Ave_sign  =  Obs_eq(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo
         

          do n =  1, Two_S
            I = L_Bath + n
            do  m =  1,Two_S
              J = L_Bath + m
                       ZXY  = GRC(I,J,1) * GR(I,J,2) +  GRC(I,J,2) * GR(I,J,1)
                       ZZ   = GRC(I,J,1) * GR(I,J,1) +  GRC(I,J,2) * GR(I,J,2)    + &
                            (GRC(I,I,2) - GRC(I,I,1))*(GRC(J,J,2) - GRC(J,J,1))
                       
                Obs_eq(1)%Obs_Latt(1,1,n,m) =  Obs_eq(1) %Obs_Latt(1,1,n,m) +  ZZ  *ZP*ZS
                Obs_eq(2) %Obs_Latt(1,1,n,m) =  Obs_eq(2) %Obs_Latt(1,1,n,m) +  ZXY *ZP*ZS
                Obs_eq(3)%Obs_Latt(1,1,n,m)=Obs_eq(3)%Obs_Latt(1,1,n,m)+ (2.d0*ZXY + ZZ)*ZP*ZS/3.d0
              
                


            enddo
          enddo 


       endif

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
          Complex (Kind=Kind(0.d0)) :: ZP, ZS, Z, ZZ, ZXY 
          Integer ::  n,m,  I,J,  I_c, I_f, J_c, J_f
          ! Add local variables as needed

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          ZS = ZS * Mc_step_weight

          ! Compute observables

         If (N_FL == 1 )  then    
          If (NT == 0 ) then
            Do I = 1,Size(Obs_tau,1)
               Obs_tau(I)%N         =  Obs_tau(I)%N + 1
               Obs_tau(I)%Ave_sign  =  Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
            Enddo
         Endif   
          do n =  1, Two_S
            I = L_Bath + n
            do  m =  1,Two_S
              J = L_Bath + m
              Z  =  - 2.d0* (G0T(J,I,1) * GT0(I,J,1))
              Obs_tau(1)%Obs_Latt(1,NT+1,n,m) =  Obs_tau(1)%Obs_Latt(1,NT+1,n,m) +   Z*ZP*ZS
            enddo
          enddo 
          Do I_c = 1,L_Bath
            do n = 1,Two_S
              I_f = L_Bath + n
              do J_c = 1,L_Bath
                do m = 1,Two_S
                  J_f = L_Bath + m
                  Z = Predefined_Obs_Cotunneling(I_c, I_f, J_c, J_f,  GT0,G0T,G00,GTT, N_SUN, N_FL) 
                  Z  =  Z  * 2.d0* g(i_c)*delta_eps(i_c) * g(j_c) * delta_eps(j_c) /real(L_Bath,kind(0.d0))
                  Obs_tau(2)%Obs_Latt(1,NT+1,1,1) =  Obs_tau(2)%Obs_Latt(1,NT+1,1,1) +   Z*ZP*ZS
                enddo
              enddo
            enddo
          enddo

    else
          If (NT == 0 ) then
            Do I = 1,Size(Obs_tau,1)
               Obs_tau(I)%N         =  Obs_tau(I)%N + 1
               Obs_tau(I)%Ave_sign  =  Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
            Enddo
         Endif

          do n =  1, Two_S
            I = L_Bath + n
            do  m =  1,Two_S
              J = L_Bath + m
              ZZ  = (( (GTT(I,I,1) -  GTT(I,I,2) ) * ( G00(J,J,1)  -  G00(J,J,2) )  -  G0T(J,I,1) * GT0(I,J,1)  -  G0T(J,I,2) * GT0(I,J,2)) )
              ZXY=  -  G0T(J,I,1) * GT0(I,J,2)  -  G0T(J,I,2) * GT0(I,J,1) 
              Obs_tau(1)%Obs_Latt(1,NT+1,n,m) =  Obs_tau(1)%Obs_Latt(1,NT+1,n,m) +   ZZ*ZP*ZS
              Obs_tau(2)%Obs_Latt(1,NT+1,n,m) =  Obs_tau(2)%Obs_Latt(1,NT+1,n,m) +   ZXY*ZP*ZS
              Obs_tau(3)%Obs_Latt(1,NT+1,n,m) =  Obs_tau(3)%Obs_Latt(1,NT+1,n,m) +   (2.d0*ZXY + ZZ)*ZP*ZS/3.d0

            enddo
          enddo 
          Do I_c = 1,L_Bath
            do n = 1,Two_S
              I_f = L_Bath + n
              do J_c = 1,L_Bath
                do m = 1,Two_S
                  J_f = L_Bath + m
                  Z = Predefined_Obs_Cotunneling(I_c, I_f, J_c, J_f,  GT0,G0T,G00,GTT, N_SUN, N_FL) 
                  Z  =  Z  * 2.d0* g(i_c)*delta_eps(i_c) * g(j_c) * delta_eps(j_c) /real(L_Bath,kind(0.d0))
                  Obs_tau(4)%Obs_Latt(1,NT+1,1,1) =  Obs_tau(4)%Obs_Latt(1,NT+1,1,1) +   Z*ZP*ZS
                enddo
              enddo
            enddo
          enddo

  endif

        end Subroutine OBSERT

!-------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief
!> Sets  the  prefactor  for  flavor  symmetry
!--------------------------------------------------------------------

       Subroutine  Set_Prefactor(Particle_Hole)

          Implicit none

          Logical, Intent(out) ::  Particle_Hole

          Integer i,  ix, iy, nc, nc1, no, no1,  i_f, i_c
          Real (Kind=Kind(0.d0))  ::  ic_p(2), X

          Allocate(Prefactor(Ndim))
          Particle_hole = .true.
          
          Prefactor  = 0
          Do  i =  1, Ndim
             Prefactor(i)  = -1
            ! if (mod(i,2)  == 0 )   Prefactor(i)  = 1
          enddo
   

end Subroutine  Set_Prefactor

!--------------------------------------------------------------------
!> @brief
!> Reconstructs dependent flavors of the configuration's weight.
!> @details
!> This has to be overloaded in the Hamiltonian submodule.
!--------------------------------------------------------------------
        subroutine weight_reconstruction(weight)
          implicit none
          complex (Kind=Kind(0.d0)), Intent(inout) :: weight(:)
          
          weight(nf_reconst) = conjg(Weight(nf_calc))  
          
        end subroutine weight_reconstruction




!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reconstructs dependent flavors of equal time Greens function
!> @details
!> This has to be overloaded in the Hamiltonian submodule.
!> @param [INOUT] Gr   Complex(:,:,:)
!> \verbatim
!>  Green function: Gr(I,J,nf) = <c_{I,nf } c^{dagger}_{J,nf } > on time slice ntau
!> \endverbatim
!-------------------------------------------------------------------
        subroutine GR_reconstruction(GR)

          Implicit none

          Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: GR(Ndim,Ndim,N_FL)
          Integer :: I,J
          complex (kind=kind(0.d0))  ::  ZZ
          Real    (kind=kind(0.d0))  :: X
        
          Do J = 1,Ndim
             Do I = 1,Ndim
                X =  real(Prefactor(I)*Prefactor(J), kind(0.d0))
                ZZ=cmplx(0.d0,0.d0,Kind(0.d0))
                if (I==J) ZZ=cmplx(1.d0,0.d0,Kind(0.d0))
                GR(I,J,nf_reconst) = ZZ - X*conjg(GR(J,I,nf_calc))
             Enddo
          Enddo
   
      end Subroutine GR_reconstruction



!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reconstructs dependent flavors of time displaced Greens function G0T and GT0
!> @details
!> This has to be overloaded in the Hamiltonian submodule.
!> @param [INOUT] GT0, G0T,  Complex(:,:,:)
!> \verbatim
!>  Green functions:
!>  GT0(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(0  )>
!>  G0T(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(tau)>
!> \endverbatim
!-------------------------------------------------------------------
      Subroutine GRT_reconstruction(GT0, G0T)
        Implicit none

        Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: GT0(Ndim,Ndim,N_FL), G0T(Ndim,Ndim,N_FL)
        Integer :: I,J
        real (kind=kind(0.d0)) :: X
        
        Do J = 1,NDIM
           Do I = 1,NDIM 
              X =  real(Prefactor(I)*Prefactor(J), kind(0.d0))
              G0T(I,J,nf_reconst) = -X*conjg(GT0(J,I,nf_calc))
              GT0(I,J,nf_reconst) = -X*conjg(G0T(J,I,nf_calc))
           enddo
         enddo
       end Subroutine GRT_reconstruction    

        
    end submodule ham_Kondo_Impurity_Wilson_smod
