      Module Hamiltonian

      Use Operator_mod
      Use Lattices_v3 
      Use MyMats 
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      Use Observables

      
      Implicit none

!>    Public variables. Have to be set by user 
      Type (Operator), dimension(:,:), allocatable  :: Op_V
      Type (Operator), dimension(:,:), allocatable  :: Op_T
      Integer, allocatable :: nsigma(:,:)
      Integer              :: Ndim,  N_FL,  N_SUN,  Ltrot
!>    Defines MPI communicator 
      Integer              :: Group_Comm


!>    Privat variables 
      Type (Lattice),       private :: Latt 
      Integer,              private :: L1, L2
      real (Kind=Kind(0.d0)),        private :: ham_T , ham_U,  Ham_chem, Ham_h, Ham_J, Ham_xi, Ham_F
      real (Kind=Kind(0.d0)),        private :: Dtau, Beta
      Character (len=64),   private :: Model, Lattice_type
      Logical,              private :: One_dimensional
      Integer,              private :: N_coord, Norb
      Integer, allocatable, private :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell
      Complex (kind=Kind(0.d0)),  private, allocatable :: ZKRON(:,:)



!>    Privat Observables
      Type (Obser_Vec ),  private, dimension(:), allocatable ::   Obs_scal
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_eq
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_tau
      
!>    Storage for the Ising action
      Real (Kind=Kind(0.d0)), private :: DW_Ising_tau(-1:1), DW_Ising_Space(-1:1), DW_Ising_Flux(-1:1,-1:1)
      Integer, allocatable  , private :: L_bond(:,:), L_bond_inv(:,:), Ising_nnlist(:,:)

    contains 


      Subroutine Ham_Set

          Implicit none

#if defined (MPI) || defined(TEMPERING)
          include 'mpif.h'
#endif   

          integer :: ierr
          Character (len=64) :: file1
          
          NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model

          NAMELIST /VAR_Hubbard/  ham_T, ham_chem, ham_U, Dtau, Beta


          NAMELIST /VAR_Hub_Ising/  ham_T, ham_chem, ham_U, Dtau, Beta, &
               &                    Ham_h, Ham_J, Ham_xi, Ham_F
#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
          !if ( irank_g == 0 )   write(6,*) "Mpi Test", igroup, isize_g
#endif
          

#ifdef MPI
          If (Irank == 0 ) then
#endif
             ! These are global variables
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
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
#endif
          Call Ham_latt

          If ( Model == "Hubbard_Mz") then
             N_FL  = 2
             N_SUN = 1
          elseif  ( Model == "Hubbard_SU2" ) then
             N_FL  = 1
             N_SUN = 2
          elseif  ( Model == "Hubbard_SU2_Ising" ) then
             N_FL  = 1
             N_SUN = 2
             If ( Lattice_type == "Honeycomb" ) then 
                Write(6,*) "Hubbard_SU2_Ising is only implemented for a square lattice"
                Stop
             Endif
          else
             Write(6,*) "Model not yet implemented!"
             Stop
          endif

#if defined(MPI) 
          If (Irank_g == 0 ) then
#endif
             File1 = "parameters"
#if defined(TEMPERING) 
             write(File1,'(A,I0,A)') "Temp_",igroup,"/parameters"
#endif
             OPEN(UNIT=5,FILE=File1,STATUS='old',ACTION='read',IOSTAT=ierr)
             If ( Model == "Hubbard_Mz"        )  READ(5,NML=VAR_Hubbard)
             If ( Model == "Hubbard_SU2"       )  READ(5,NML=VAR_Hubbard)
             If ( Model == "Hubbard_SU2_Ising" )  READ(5,NML=VAR_Hub_Ising)
             CLOSE(5)
#ifdef MPI
          endif
          CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_chem ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_U    ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(Dtau     ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(Beta     ,1,MPI_REAL8,0,Group_Comm,ierr)
          If ( Model == "Hubbard_SU2_Ising" ) then
             CALL MPI_BCAST(Ham_xi   ,1,MPI_REAL8,0,Group_Comm,ierr)
             CALL MPI_BCAST(Ham_J    ,1,MPI_REAL8,0,Group_Comm,ierr)
             CALL MPI_BCAST(Ham_h    ,1,MPI_REAL8,0,Group_Comm,ierr)
             CALL MPI_BCAST(Ham_F    ,1,MPI_REAL8,0,Group_Comm,ierr)
          Endif
#endif
           
           Call Ham_hop
           Ltrot = nint(beta/dtau)
           
           If  ( Model == "Hubbard_SU2_Ising" )  Call Setup_Ising_action

#if defined(TEMPERING)
           write(File1,'(A,I0,A)') "Temp_",igroup,"/info"
#else
           write(File1,'(A,I0)') "info"
#endif
           
#if defined(MPI) 
           If (Irank_g == 0)  then
#endif
             
              Open (Unit = 50,file=file1,status="unknown",position="append")
              Write(50,*) '====================================='
              Write(50,*) 'Model is      : ', Model 
              Write(50,*) 'Lattice is    : ', Lattice_type
              Write(50,*) '# of orbitals : ', Ndim
              Write(50,*) 'Beta          : ', Beta
              Write(50,*) 'dtau,Ltrot    : ', dtau,Ltrot
              Write(50,*) 'N_SUN         : ', N_SUN
              Write(50,*) 'N_FL          : ', N_FL
              Write(50,*) 't             : ', Ham_T
              Write(50,*) 'Ham_U         : ', Ham_U
              Write(50,*) 'Ham_chem      : ', Ham_chem
              If ( Model == "Hubbard_SU2_Ising" ) then
                 Write(50,*) 'Ham_xi        : ', Ham_xi
                 Write(50,*) 'Ham_J         : ', Ham_J
                 Write(50,*) 'Ham_h         : ', Ham_h
                 Write(50,*) 'Ham_F         : ', Ham_F
              Endif
              close(50)
#if defined(MPI) 
           endif
           
#endif
           call Ham_V
           
         end Subroutine Ham_Set
!=============================================================================
        Subroutine Ham_Latt
          Implicit none
          !Set the lattice
          
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          Integer :: I, nc, no

          If ( Lattice_type =="Square" ) then
             Norb = 1
             N_coord   = 2
             One_dimensional = .false.
             a1_p(1) =  1.0  ; a1_p(2) =  0.d0
             a2_p(1) =  0.0  ; a2_p(2) =  1.d0
             L1_p    =  dble(L1)*a1_p
             L2_p    =  dble(L2)*a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
             If ( L1 == 1 .or. L2 == 1 ) then 
                One_dimensional = .true.
                N_coord   = 1
                If (L1 == 1 ) then 
                   Write(6,*) ' For one dimensional systems set  L2 = 1 ' 
                   Stop
                endif
             endif
          elseif ( Lattice_type=="Honeycomb" ) then
             Norb    = 2
             N_coord = 3
             a1_p(1) =  1.D0   ; a1_p(2) =  0.d0
             a2_p(1) =  0.5D0  ; a2_p(2) =  sqrt(3.D0)/2.D0

             !del_p   =  (a2_p - 0.5*a1_p ) * 2.0/3.0
             L1_p    =  dble(L1) * a1_p
             L2_p    =  dble(L2) * a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )

          else
             Write(6,*) "Lattice not yet implemented!"
             Stop
          endif


          ! This is for the orbital structure.
          Ndim = Latt%N*Norb
          Allocate (List(Ndim,2), Invlist(Latt%N,Norb))
          nc = 0
          Do I = 1,Latt%N
             Do no = 1,Norb
                ! For the Honeycomb lattice no = 1,2 corresponds to the A,and B sublattices.
                nc = nc + 1
                List(nc,1) = I
                List(nc,2) = no
                Invlist(I,no) = nc 
             Enddo
          Enddo

        end Subroutine Ham_Latt

!===================================================================================           
        Subroutine Ham_hop
          Implicit none

          !Setup the hopping
          !Per flavor, the  hopping is given by: 
          !  e^{-dtau H_t}  =    Prod_{n=1}^{Ncheck} e^{-dtau_n H_{n,t}}

          Integer :: I, I1, J1, I2, n, Ncheck,nc, nc1, no

          Ncheck = 1
          allocate(Op_T(Ncheck,N_FL))
          do n = 1,N_FL
             Do nc = 1,Ncheck
                Call Op_make(Op_T(nc,n),Ndim)
                If ( Lattice_type =="Square" ) then
                   If (One_dimensional ) then 
                      DO I = 1, Latt%N
                         I1 = Latt%nnlist(I, 1, 0)
                         Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                      ENDDO
                   else
                      DO I = 1, Latt%N
                         I1 = Latt%nnlist(I,1,0)
                         I2 = Latt%nnlist(I,0,1)
                         Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I,I2) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I2,I) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                      ENDDO
                   endif
                elseif ( Lattice_type=="Honeycomb" ) then
                   DO I = 1, Latt%N
                      do no = 1,Norb
                         I1 = Invlist(I,no)
                         Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                      enddo
                      I1 = Invlist(I,1)
                      Do nc1 = 1,N_coord
                         select case (nc1)
                         case (1)
                            J1 = invlist(I,2) 
                         case (2)
                            J1 = invlist(Latt%nnlist(I,1,-1),2) 
                         case (3)
                            J1 = invlist(Latt%nnlist(I,0,-1),2) 
                         case default
                            Write(6,*) ' Error in  Ham_Hop '  
                         end select
                         Op_T(nc,n)%O(I1,J1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                         Op_T(nc,n)%O(J1,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                      Enddo
                   Enddo
                endif

                Do I = 1,Ndim
                   Op_T(nc,n)%P(i) = i 
                Enddo
                if ( abs(Ham_T) < 1.E-6 .and.  abs(Ham_chem) < 1.E-6 ) then 
                   Op_T(nc,n)%g = 0.d0
                else
                   Op_T(nc,n)%g = -Dtau
                endif
                Op_T(nc,n)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
                Call Op_set(Op_T(nc,n)) 
             enddo
          enddo
        end Subroutine Ham_hop

!===================================================================================           
        Subroutine Ham_V
          
          Implicit none 
          
          Integer :: nf, I, I1, I2,  nc, nc1,  J
          Real (Kind=Kind(0.d0)) :: X

          
          If (Model == "Hubbard_SU2")  then
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
          Elseif (Model == "Hubbard_Mz")  then
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
          Elseif  ( Model == "Hubbard_SU2_Ising" ) then
             If ( Abs(Ham_U) > 1.d-6) then
                Allocate(Op_V(3*Ndim,N_FL))
             else
                Allocate(Op_V(2*Ndim,N_FL))
             endif
             do nf = 1,N_FL
                do i  =  1, N_coord*Ndim
                   call Op_make(Op_V(i,nf),2)
                enddo
                If ( Abs(Ham_U) > 1.d-6) then
                   do i  = N_coord*Ndim +1 ,  N_coord*Ndim + Ndim ! For Hubbard
                      Call Op_make(Op_V(i,nf),1)
                   enddo
                endif
             enddo
             Do nc = 1,Ndim*N_coord   ! Runs over bonds.  Coordination number = 2.
                                      ! For the square lattice Ndim = Latt%N
                I1 = L_bond_inv(nc,1)
                if ( L_bond_inv(nc,2)  == 1 ) I2 = Latt%nnlist(I1,1,0)
                if ( L_bond_inv(nc,2)  == 2 ) I2 = Latt%nnlist(I1,0,1)
                Op_V(nc,1)%P(1) = I1
                Op_V(nc,1)%P(2) = I2
                Op_V(nc,1)%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.D0)) 
                Op_V(nc,1)%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.D0)) 
                Op_V(nc,1)%g      = cmplx(-dtau*Ham_xi,0.D0,kind(0.D0))
                Op_V(nc,1)%alpha  = cmplx(0d0,0.d0, kind(0.D0)) 
                Op_V(nc,1)%type   = 1
                Call Op_set( Op_V(nc,1) )
             Enddo

             If ( Abs(Ham_U) > 1.d-6) then
                Do i = 1,Ndim
                   nc1 = N_coord*Ndim + i
                   Op_V(nc1,1)%P(1)   = i
                   Op_V(nc1,1)%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
                   Op_V(nc1,1)%g      = sqrt(cmplx(-dtau*ham_U/(DBLE(N_SUN)), 0.D0, kind(0.D0)))
                   Op_V(nc1,1)%alpha  = cmplx(-0.5d0,0.d0, kind(0.d0))
                   Op_V(nc1,1)%type   = 2
                   Call Op_set( Op_V(nc1,1) )
                Enddo
             Endif
          Endif
        end Subroutine Ham_V

!===================================================================================           
        Real (Kind=Kind(0.d0)) function S0(n,nt)  
          Implicit none
          Integer, Intent(IN) :: n,nt
          Integer :: nt1,I, F1,F2,I1,I2,I3
          
          !> Ratio for local spin-flip
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

             ! Magnetic flux term
             I1 = L_bond_inv(n,1)
             if ( L_bond_inv(n,2) == 1 ) then
                !     I2
                !     I1 I3
                I2 = Latt%nnlist(I1,0,1 )
                I3 = Latt%nnlist(I1,1,0 )
                F1 = nsigma(n,nt)*nsigma(L_bond(I1,2),nt)* nsigma(L_bond(I2,1),nt)*nsigma(L_bond(I3,2),nt)
                !     I1  
                !     I2 I3 
                I2 = Latt%nnlist(I1,0,-1)
                I3 = Latt%nnlist(I1,1,-1)
                F2 = nsigma(n,nt)*nsigma(L_bond(I2,1),nt)* nsigma(L_bond(I2,2),nt)*nsigma(L_bond(I3,2),nt)
             else
                !    I3
                !    I2  I1
                I2 = Latt%nnlist(I1,-1,0 )
                I3 = Latt%nnlist(I1,-1,1 )
                F1 = nsigma(n,nt)*nsigma(L_bond(I2,1),nt)* nsigma(L_bond(I2,2),nt)*nsigma(L_bond(I3,1),nt)
                !    I2
                !    I1  I3
                I2 = Latt%nnlist(I1,0,1)
                I3 = Latt%nnlist(I1,1,0)
                F2 = nsigma(n,nt)*nsigma(L_bond(I1,1),nt)* nsigma(L_bond(I2,1),nt)*nsigma(L_bond(I3,2),nt)
             endif
             S0 = S0*DW_Ising_Flux(F1,F2)
          endif
          
        end function S0

!===================================================================================           
        Subroutine Global_move(T0_Proposal_ratio,nsigma_old,size_clust)
          !>  The input is the field nsigma declared in this module. This routine generates a 
          !>  global update with  and returns the propability  
          !>  T0_Proposal_ratio  =  T0( sigma_out-> sigma_in ) /  T0( sigma_in -> sigma_out)  
          !>   
          
          Implicit none
          Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio, size_clust
          Integer, dimension(:,:),  allocatable, intent(in)  :: nsigma_old
          !> nsigma_old contains a copy of nsigma upon entry
          
          !> Local
          Integer :: I,I1, nt,nt1,nt_end,  n, ns_old, n1,n2,n3,n4, dx, dy,dt,sx,sy,st, nc, nb, n_kinks
          Integer :: ndis, ndir,m, n_L, n_R, nw, nw_r,N_samples
          Real (Kind=Kind(0.d0)) :: Weight, Ratio, X, Y
          Integer, allocatable ::  Kinks(:), Pos(:)
          Real (Kind=Kind(0.d0)), allocatable :: W(:), W1(:)

          Logical :: Test

          Test = .false.
          If (Model == "Hubbard_SU2_Ising" ) then
             if ( ranf_wrap() > 0.5d0 ) then
                ! Move kinks along the imaginary time axis.
                if (Test) write(6,*) '********Global *******'
                Allocate ( Kinks(Ltrot), W(0:Ltrot), Pos(0:Ltrot), W1(0:Ltrot) )
                Kinks   = 0
                Pos     = 0
                W       = 0.d0
                n_kinks = 0
                Kinks   = 0
                I = nranf(Latt%N)
                nb= nranf(N_coord)
                do nt = 1,Ltrot
                   nt1 = nt + 1
                   if ( nt == Ltrot ) nt1 = 1
                   if ( nsigma(L_bond(I,nb),nt)*nsigma(L_bond(I,nb),nt1) == -1 ) then
                      n_kinks   = n_kinks + 1
                      Kinks(n_kinks) = nt
                   endif
                enddo
                
                If (test) then
                   Write(6,*) "----:", n_kinks
                   do n = 1,n_kinks
                      Write(6,*) Kinks(n)
                   enddo
                endif
                
                if ( n_kinks > 1 ) then
                   n = nranf(n_kinks) ! The kink
                   Pos(0) = Kinks(n)
                   W  (0) = 1.d0
                   
                   if (test) WRITE(6,*) "The Kink:",N
                   nc = 0
                   n_L = n -1   ! The left kink
                   if (n_L == 0) n_L = n_kinks
                   n_R = n + 1  ! The right kink
                   if (n_R > n_kinks) n_R = n_R - n_kinks
                   nt = Kinks(n)
                   X  = 1.d0
                   nw = 0
                   do
                      nt  = nt + 1
                      if (nt > Ltrot ) nt = 1
                      if ( nt == Kinks(n_R) ) exit
                      X = X*S0(L_bond(I,nb),nt)
                      nw = nw + 1
                      !Write(6,*) "R,nw: ", nw, S0(L_bond(I,nb),nt)
                      W  (nw) = X
                      Pos(nw) = nt
                      nsigma(L_bond(I,nb),nt) = -  nsigma(L_bond(I,nb),nt)
                   enddo
                   nsigma =  nsigma_old
                   nw_r = nw
                   nt = Kinks(n) 
                   X  = 1.d0
                   do
                      nt = nt - 1
                      if ( nt == 0  )   nt = Ltrot
                      if ( nt == Kinks(n_L) ) exit
                      nt1 = nt+1
                      if (nt1 > Ltrot) nt1 = 1
                      X = X*S0(L_bond(I,nb),nt1)
                      nw = nw + 1
                      !Write(6,*) "L,nw: ", nw, S0(L_bond(I,nb),nt1)
                      W(nw)   = X
                      Pos(nw) = nt
                      nsigma(L_bond(I,nb),nt1) = -  nsigma(L_bond(I,nb),nt1)
                   enddo
                   nsigma =  nsigma_old
                   
                   ! Normalize
                   X = 0.d0
                   do  m = 0,Nw
                      X = X + W(m)
                   enddo
                   W = W/X
                   ! Sample
                   X = ranf_wrap()
                   Y = 0.d0
                   do m = 0,Nw
                      Y = Y + W(m)
                      if (X < Y) exit
                   enddo
                   nt_end = Pos(m)
                   
                   if (test) Write(6,*) 'Move: ', Pos(0), " to ", Pos(m)
                   if ( m == 0 ) then
                      T0_Proposal_ratio = 0.d0
                   elseif ( m <= nw_r )  then
                      nt = Kinks(n)
                      if (test) Write(6,*) 'Up'
                      do
                         if ( nt == nt_end ) exit
                         nt  = nt + 1
                         if (nt > Ltrot ) nt = 1
                         nsigma(L_bond(I,nb),nt) = -  nsigma(L_bond(I,nb),nt)
                      enddo
                      T0_Proposal_ratio = W(0)/W(m)
                   else
                      nt = Kinks(n) 
                      if (test) Write(6,*) 'do'
                      do
                         if ( nt ==  nt_end ) exit
                         nsigma(L_bond(I,nb),nt) = -  nsigma(L_bond(I,nb),nt)
                         nt = nt - 1
                         if ( nt == 0  )   nt = Ltrot
                      enddo
                      T0_Proposal_ratio = W(0)/W(m)
                   end if
                   if (test)  then
                      ! Check
                      n_kinks = 0
                      do nt = 1,Ltrot
                         nt1 = nt + 1
                         if ( nt == Ltrot ) nt1 = 1
                         if ( nsigma(L_bond(I,nb),nt)*nsigma(L_bond(I,nb),nt1) == -1 ) then
                            n_kinks   = n_kinks + 1
                            Kinks(n_kinks) = nt
                         endif
                      enddo
                      Write(6,*) "----:", n_kinks
                      do n = 1,n_kinks
                         Write(6,*) Kinks(n)
                      enddo
                   endif
                   deallocate (Kinks, W, Pos,W1)
                Endif
             else
                I = nranf(Latt%N)
                n1  = L_bond(I,1)
                n2  = L_bond(I,2)
                n3  = L_bond(Latt%nnlist(I,-1,0),1)
                n4  = L_bond(Latt%nnlist(I,0,-1),2)
                do nt = 1,Ltrot
                   nsigma(n1,nt) = -nsigma(n1,nt)
                   nsigma(n2,nt) = -nsigma(n2,nt)
                   nsigma(n3,nt) = -nsigma(n3,nt)
                   nsigma(n4,nt) = -nsigma(n4,nt)
                enddo
                Ratio = Delta_S0_global(Nsigma_old)
                Weight = 1.d0 - 1.d0/(1.d0+Ratio)
                If ( Weight < ranf_wrap() ) Then
                   T0_Proposal_ratio = 0.d0
                   nsigma            = nsigma_old
                else
                   T0_Proposal_ratio = 1.d0/Ratio
                endif
             endif
          endif

        End Subroutine Global_move
!===================================================================================           
        Real (Kind=kind(0.d0)) Function Delta_S0_global(Nsigma_old)

          !>  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none 

          !> Arguments
          Integer, dimension(:,:), allocatable, intent(IN) :: Nsigma_old
          !> Local 
          Integer :: I,n,n1,n2,n3,n4,nt,nt1, nc_F, nc_J, nc_h_p, nc_h_m
         

          Delta_S0_global = 1.d0
          If ( Model == "Hubbard_SU2_Ising" ) then
             nc_F = 0
             nc_J = 0
             nc_h_p = 0
             nc_h_m = 0
             Do I = 1,Latt%N
                n1  = L_bond(I,1)
                n2  = L_bond(Latt%nnlist(I,1,0),2)
                n3  = L_bond(Latt%nnlist(I,0,1),1)
                n4  = L_bond(I,2)
                do nt = 1,Ltrot
                   nt1 = nt +1 
                   if (nt == Ltrot) nt1 = 1
                   if (nsigma(n1,nt) == nsigma(n1,nt1) ) then 
                      nc_h_p = nc_h_p + 1
                   else
                      nc_h_m = nc_h_m + 1
                   endif
                   if (nsigma_old(n1,nt) == nsigma_old(n1,nt1) ) then 
                      nc_h_p = nc_h_p - 1
                   else
                      nc_h_m = nc_h_m - 1
                   endif

                   if (nsigma(n4,nt) == nsigma(n4,nt1) ) then 
                      nc_h_p = nc_h_p + 1
                   else
                      nc_h_m = nc_h_m + 1
                   endif
                   if (nsigma_old(n4,nt) == nsigma_old(n4,nt1) ) then 
                      nc_h_p = nc_h_p - 1
                   else
                      nc_h_m = nc_h_m - 1
                   endif
                   
                   nc_F = nc_F + nsigma    (n1,nt)*nsigma    (n2,nt)*nsigma    (n3,nt)*nsigma    (n4,nt)  &
                        &      - nsigma_old(n1,nt)*nsigma_old(n2,nt)*nsigma_old(n3,nt)*nsigma_old(n4,nt) 
                   
                   nc_J = nc_J + nsigma(n1,nt)*nsigma(n2,nt) + &
                        &        nsigma(n2,nt)*nsigma(n3,nt) + &
                        &        nsigma(n3,nt)*nsigma(n4,nt) + &
                        &        nsigma(n4,nt)*nsigma(n1,nt) - &
                        &        nsigma_old(n1,nt)*nsigma_old(n2,nt) - &
                        &        nsigma_old(n2,nt)*nsigma_old(n3,nt) - &
                        &        nsigma_old(n3,nt)*nsigma_old(n4,nt) - &
                        &        nsigma_old(n4,nt)*nsigma_old(n1,nt) 
                   
                enddo
             enddo
             Delta_S0_global = ( sinh(Dtau*Ham_h)**nc_h_m ) * (cosh(Dtau*Ham_h)**nc_h_p) * &
                  &            exp( -Dtau*(Ham_F*real(nc_F,kind(0.d0)) -  Ham_J*real(nc_J,kind(0.d0))))
          endif
        end Function Delta_S0_global
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
          DO n = -1,1,2
             do n1 = -1,1,2
                DW_Ising_Flux(n,n1) = exp( Dtau*Ham_F*(dble(n) +  dble(n1) ))/ exp(  -Dtau*Ham_F*(dble(n) +  dble(n1) ))
             enddo
          enddo
          
        End Subroutine Setup_Ising_action
!===================================================================================           
        Subroutine  Alloc_obs(Ltau) 

          Implicit none
          Integer, Intent(In) :: Ltau
          Integer    ::  i, N, Ns,Nt,No
          Character (len=64) ::  Filename

          
          Allocate (ZKRON(Ndim,Ndim))
          ZKRON = cmplx(0.d0,0.d0,Kind(0.d0))
          Do I = 1,Ndim
             ZKRON(I,I) = cmplx(1.d0,0.d0,kind(0.d0))
          Enddo
          

          ! Scalar observables
          Allocate ( Obs_scal(14) )
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
             case (5)
                N = 1;   Filename ="Flux"
             case (6)
                N = 1;   Filename ="X"
             case (7)
                N = L1/2       
                If (L2 > L1 ) N = L2/2; Filename ="Wilson"
             case (8)
                N = L1/2 + L2/2;        Filename ="Vison"
             case (9)
                N = L1/2       
                If (L2 > L1 ) N = L2/2; Filename ="Prop"
             case (10)
                N = 1      
                Filename ="Q"
             case (11)
                N = 1      
                Filename ="Star"
             case (12)
                N = 1 
                Filename ="QT"
             case (13)
                N = 1 
                Filename ="QTT"
             case (14)
                N = 1 
                Filename ="QStag"
             case default
                Write(6,*) ' Error in Alloc_obs '  
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
          enddo


          ! Equal time correlators
          Allocate ( Obs_eq(7) )
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
                Ns = Latt%N;  No = Norb;  Filename ="Flux"
             case (6)
                Ns = Latt%N;  No = N_coord;  Filename ="Kin"
             case (7)
                Ns = Latt%N;  No = N_coord;  Filename ="Dimer"
             case default
                Write(6,*) ' Error in Alloc_obs '  
             end select
             Nt = 1
             Call Obser_Latt_make(Obs_eq(I),Ns,Nt,No,Filename)
          enddo

          If (Ltau == 1) then 
             ! Equal time correlators
             Allocate ( Obs_tau(7) )
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
                case (5)
                   Ns = Latt%N;  No = N_coord;  Filename ="X"
                case (6)
                   Ns = Latt%N;  No = N_coord;  Filename ="Kin"
                case (7)
                   Ns = Latt%N;  No = N_coord;  Filename ="Dimer"
                case default
                   Write(6,*) ' Error in Alloc_obs '  
                end select
                Nt = Ltrot+1
                Call Obser_Latt_make(Obs_tau(I),Ns,Nt,No,Filename)
             enddo
          endif
        end Subroutine Alloc_obs

!========================================================================
        Subroutine Obser(GR,Phase,Ntau)
          
          Implicit none

          
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          
          !Local 
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK, G(4,4)
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS, ZN
          Complex (Kind=Kind(0.d0)) :: ZQ, ZSTAR, ZQT, ZQTT, ZQP
          Integer :: I,J, imj, nf, dec, I1, I2,I3,I4, J1, no_I, no_J,  iFlux_tot,  &
               &     no, no1, ntau1, L_Vison, L_Wilson, n, nx,ny
          Real (Kind=Kind(0.d0)) :: X_ave, X, XI1,XI2,XI3,XI4, Y


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
          ! GRC(i,j,nf) = < c^{dagger}_{i,nf } c_{j,nf } >

          ! Compute scalar observables. 
          Do I = 1,Size(Obs_scal,1)
             Obs_scal(I)%N         =  Obs_scal(I)%N + 1
             Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo
             

          Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do J = 1,Ndim
                Zkin = Zkin + sum(Op_T(1,nf)%O(:, j)*Grc(:, j, nf))
             ENddo
          Enddo
          Zkin = Zkin * dble(N_SUN)
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *ZP* ZS


          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising" ) then
            dec = 1
          else
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


          iFlux_tot = 0
          Do I = 1, Ndim
             iFlux_tot = iFlux_tot + iFlux(I,Ntau)
          Enddo
          Obs_scal(5)%Obs_vec(1)  =   Obs_scal(5)%Obs_vec(1) + cmplx(dble(iFlux_tot),0.d0,kind(0.d0))*ZP*ZS


          X_ave = 0.d0
          ntau1 = ntau + 1
          If (ntau == Ltrot)  ntau1 = 1
          Do I = 1,Latt%N
             do no = 1,2
                X_ave = X_ave + DW_Ising_tau( nsigma(L_bond(I,no),ntau)*nsigma(L_bond(I,no),ntau1) )
             Enddo
          Enddo
          Obs_scal(6)%Obs_vec(1) = Obs_scal(6)%Obs_vec(1) + cmplx(X_ave,0.d0,kind(0.d0))*ZP*ZS

          Do I = 1,Latt%N
             do L_Wilson = 1,Size(Obs_scal(7)%Obs_vec,1)
                X  = 1.d0
                I1 = I
                I2 = I
                Do n = 1,L_Wilson
                   X = X *nsigma(L_bond(I1,1),ntau)
                   I1 = Latt%nnlist(I1,1,0)
                   X = X *nsigma(L_bond(I2,2),ntau)
                   I2 = Latt%nnlist(I2,0,1)
                enddo
                Do n = 1,L_Wilson
                   X = X *nsigma(L_bond(I1,2),ntau)
                   I1 = Latt%nnlist(I1,0,1)
                   X = X* nsigma(L_bond(I2,1),ntau)
                   I2 = Latt%nnlist(I2,1,0)
                enddo
                Obs_scal(7)%Obs_vec(L_Wilson) = Obs_scal(7)%Obs_vec(L_Wilson) + cmplx(X,0.d0,kind(0.d0)) *  ZP * ZS

                X  = 1.d0
                I1 = I
                Do n = 1,L_Wilson/2
                   X = X *nsigma(L_bond(I1,2),ntau)
                   I1 = Latt%nnlist(I1,0,1)
                enddo
                Do n = 1,L_Wilson
                   X = X *nsigma(L_bond(I1,1),ntau)
                   I1 = Latt%nnlist(I1,1,0)
                enddo
                Do n = 1,L_Wilson/2
                   I1 = Latt%nnlist(I1,0,-1)
                   X = X *nsigma(L_bond(I1,2),ntau)
                enddo
                Z = cmplx(0.d0,0.d0,kind(0.d0))
                Do nf = 1,N_FL
                   Z = Z  + cmplx(X,0.d0,kind(0.d0))* Grc(I, I1, nf)
                enddo
                Obs_scal(9)%Obs_vec(L_Wilson)   = Obs_scal(9)%Obs_vec(L_Wilson) +  Z *  ZP * ZS
             Enddo
          Enddo

          
          Do I = 1, Latt%N
             I1 = I
             X  = 1.d0
             DO L_Vison = 1,L1/2
                X = X*DW_Ising_tau( nsigma(L_bond(I1 ,2),ntau)*nsigma(L_bond(I1 ,2),ntau1) )
                Obs_scal(8)%Obs_vec(L_Vison) =   Obs_scal(8)%Obs_vec(L_Vison) + cmplx(X,0.d0,kind(0.d0)) * ZP * ZS
                I1 = Latt%nnlist(I1,1,0)
             Enddo
             I1 = Latt%nnlist(I1,0,1)
             DO L_Vison = L1/2 + 1, L1/2 + L2/2
                X = X*DW_Ising_tau( nsigma(L_bond(I1 ,1),ntau)*nsigma(L_bond(I1 ,1),ntau1) )
                Obs_scal(8)%Obs_vec(L_Vison) =   Obs_scal(8)%Obs_vec(L_Vison) + cmplx(X,0.d0,kind(0.d0)) * ZP * ZS
                I1 = Latt%nnlist(I1,0,1)
             Enddo
           Enddo

           ZQ    = cmplx(0.d0,0.d0,kind(0.d0))
           ZQP   = cmplx(0.d0,0.d0,kind(0.d0))
           ZQT   = cmplx(0.d0,0.d0,kind(0.d0))
           ZQTT  = cmplx(0.d0,0.d0,kind(0.d0))
           ZStar = cmplx(0.d0,0.d0,kind(0.d0))
           If (ntau == Ltrot)  ntau1 = 1
           Do I = 1,Latt%N
              If ( mod(N_SUN,2) == 0 ) then
                 X = 1.d0
                 Y = (-1.d0)**( Latt%List(I,1) +  Latt%List(I,2) )
              else
                 X = (-1.d0)**( Latt%List(I,1) +  Latt%List(I,2) )
              endif
              I1 = Latt%nnlist(I, 1, 0)
              I2 = Latt%nnlist(I, 0, 1)
              I3 = Latt%nnlist(I,-1, 0)
              I4 = Latt%nnlist(I, 0,-1)
              XI1 = DW_Ising_tau( nsigma(L_bond(I ,1),ntau)*nsigma(L_bond(I ,1),ntau1) )
              XI2 = DW_Ising_tau( nsigma(L_bond(I ,2),ntau)*nsigma(L_bond(I ,2),ntau1) )
              XI3 = DW_Ising_tau( nsigma(L_bond(I3,1),ntau)*nsigma(L_bond(I3,1),ntau1) )
              XI4 = DW_Ising_tau( nsigma(L_bond(I4,2),ntau)*nsigma(L_bond(I4,2),ntau1) )
              Z    =  ( cmplx(1.d0,0.d0,kind(0.d0)) - cmplx(2.d0,0.d0,kind(0.d0))*GRC(I,I,1) )**(N_SUN)
              ZQ   = ZQ    + cmplx(X*XI1*XI2*XI3*XI4,0.d0,kind(0.d0)) * Z
              ZQP  = ZQP   + cmplx(Y*X*XI1*XI2*XI3*XI4,0.d0,kind(0.d0)) * Z
              ZQT  = ZQT   + cmplx(X*XI1*XI2*XI3*XI4,0.d0,kind(0.d0))  +  Z
              ZQTT = ZQTT  + cmplx(X*XI1*XI2,0.d0,kind(0.d0))  +  cmplx(X*XI3*XI4,0.d0,kind(0.d0))*Z
              ZStar= ZStar + cmplx(XI1*XI2*XI3*XI4  ,0.d0,kind(0.d0))
           Enddo
           Obs_scal(10)%Obs_vec(1) = Obs_scal(10)%Obs_vec(1) + ZQ    * ZP*ZS
           Obs_scal(11)%Obs_vec(1) = Obs_scal(11)%Obs_vec(1) + ZStar * ZP*ZS
           Obs_scal(12)%Obs_vec(1) = Obs_scal(12)%Obs_vec(1) + ZQT   *ZP*ZS
           Obs_scal(13)%Obs_vec(1) = Obs_scal(13)%Obs_vec(1) + ZQTT  *ZP*ZS
           Obs_scal(14)%Obs_vec(1) = Obs_scal(14)%Obs_vec(1) + ZQP   *ZP*ZS
           
           
           ! Compute spin-spin, Green, and den-den correlation functions  !  This is general N_SUN, and  N_FL = 1
           DO I = 1,Size(Obs_eq,1)
              Obs_eq(I)%N        = Obs_eq(I)%N + 1
              Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
           ENDDO
           If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising"  ) then 
              ZN =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
              Do I1 = 1,Ndim
                 I    = List(I1,1)
                 no_I = List(I1,2)
                 Do J1 = 1,Ndim
                    J    = List(J1,1)
                    no_J = List(J1,2)
                    imj = latt%imj(I,J)
                    ! Green
                    Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + &
                         &               ZN * GRC(I1,J1,1) *  ZP*ZS 
                    ! SpinZ
                    Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + &
                         &               ZN * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                    ! SpinXY
                    Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + &
                         &               ZN * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                    ! Den
                    Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J)  +  &
                         &     (    GRC(I1,I1,1) * GRC(J1,J1,1) * ZN     + &
                         &          GRC(I1,J1,1) * GR(I1,J1,1 )           &
                         &                                    ) * ZN* ZP*ZS
                    !Flux
                    Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) + &
                         &         cmplx(real(iFlux(I1,Ntau)*iFlux(J1,Ntau),kind(0.d0)),0.d0,kind(0.d0))*ZP*ZS
                 ENDDO
                 Obs_eq(4)%Obs_Latt0(no_I) =  Obs_eq(4)%Obs_Latt0(no_I) +  Z * GRC(I1,I1,1) * ZP * ZS
                 Obs_eq(5)%Obs_Latt0(no_I) =  Obs_eq(5)%Obs_Latt0(no_I) +  &
                      &   cmplx(real(iFlux(I1,Ntau),kind(0.d0)),0.d0,kind(0.d0)) * ZP * ZS
              ENDDO
              ! Dimer-Dimer, Kin-Kin correlations
              Do I = 1,Latt%N
                 do no = 1,N_coord
                    select  case(no)
                    case(1) 
                       I1 = Latt%nnlist(I,1,0)
                    case(2)
                       I1 = Latt%nnlist(I,0,1)
                    end select
                    Do J = 1,Latt%N
                       imj = latt%imj(I,J)
                       do no1 = 1,N_coord
                          if ( no1 == 1) J1 = Latt%nnlist(J,1,0)
                          if ( no1 == 2) J1 = Latt%nnlist(J,0,1)
                          G(1,1)=  GRC(I,I,1); G(1,2)= GRC(I,I1,1) ;G(1,3)= GRC(I,J,1) ;G(1,4)=GRC(I,J1,1)
                          G(2,1)= -GR(I,I1,1); G(2,2)= GRC(I1,I1,1);G(2,3)= GRC(I1,J,1);G(2,4)=GRC(I1,J1,1)
                          G(3,1)= -GR(I,J,1) ; G(3,2)=-GR(I1,J,1)  ;G(3,3)= GRC(J,J,1) ;G(3,4)=GRC(J,J1,1)
                          G(4,1)= -GR(I,J1,1); G(4,2)=-GR(I1,J1,1) ;G(4,3)=-GR(J,J1,1) ;G(4,4)=GRC(J1,J1,1)
                          If ( N_SUN == 4 ) then
#include "Bid_N4.f90"
                          endif
                          If ( N_SUN == 3 ) then
#include "Bid_N3.f90"
                          endif
                          If ( N_SUN == 2 ) then
#include "Bid_N2.f90"
                          endif
                          Obs_eq(7)%Obs_Latt(imj,1,no,no1)  = Obs_eq(7)%Obs_Latt(imj,1,no,no1)   +  Z * ZP*ZS

                          Z =  (  (GRC(I,I1,1) +  GRC(I1,I,1)) * (GRC(J,J1,1)  +  GRC(J1,J,1)) * ZN  + &
                               &   GRC(I ,J1,1)*GR(I1,J ,1) + &
                               &   GRC(I ,J ,1)*GR(I1,J1,1) + &
                               &   GRC(I1,J1,1)*GR(I ,J ,1) + &
                               &   GRC(I1,J ,1)*GR(I ,J1,1)   )* ZN * &
                               &   real(nsigma(L_bond(I ,no),ntau)*nsigma(L_bond(J ,no1),ntau), Kind(0.d0))
                          Obs_eq(6)%Obs_Latt(imj,1,no,no1)  = Obs_eq(6)%Obs_Latt(imj,1,no,no1)   +  Z * ZP*ZS
                       enddo

                    enddo
                    If ( N_SUN == 2 ) then
#include "Bid_N2_0.f90"
                    endif
                    Obs_eq(7)%Obs_Latt0(no) =  Obs_eq(7)%Obs_Latt0(no) +  Z * ZP*ZS
                    Z = (GRC(I,I1,1) +  GRC(I1,I,1))*ZN*Real(nsigma(L_bond(I ,no),ntau),kind(0.d0))
                    Obs_eq(6)%Obs_Latt0(no) =  Obs_eq(6)%Obs_Latt0(no) +  Z * ZP*ZS

                 enddo
              enddo



              ! Compute Vison-Vison correlations
!!$             Do I = 1,Latt%N
!!$                X = 1.d0
!!$                J = I
!!$                Do ny = 1,L2/2
!!$                   J   = Latt%nnlist(J,1,0)
!!$                   Do nx = 1,L1-1
!!$                      X   = X*DW_Ising_tau( nsigma(L_bond(J ,2),ntau)*nsigma(L_bond(J ,2),ntau1) )
!!$                      imj = Latt%imj(i,j)
!!$                      Obs_eq(6)%Obs_Latt(imj,1,1,1) = Obs_eq(6)%Obs_Latt(imj,1,1,1) + &
!!$                           & cmplx(X,0.d0,kind(0.d0))* ZP * ZS
!!$                      J   = Latt%nnlist(J,1,0)
!!$                   Enddo
!!$                   ! J = I + L1*a_1
!!$                   J = Latt%nnlist(J,-1,1)
!!$                   X   = X*DW_Ising_tau( nsigma(L_bond(J ,1),ntau)*nsigma(L_bond(J ,1),ntau1) )
!!$                   imj = Latt%imj(i,j)
!!$                   Obs_eq(6)%Obs_Latt(imj,1,1,1) = Obs_eq(6)%Obs_Latt(imj,1,1,1) + &
!!$                        & cmplx(X,0.d0,kind(0.d0))* ZP * ZS
!!$                   Do nx = 1,L1-1
!!$                      J   = Latt%nnlist(J,-1,0)
!!$                      X   = X*DW_Ising_tau( nsigma(L_bond(J ,2),ntau)*nsigma(L_bond(J ,2),ntau1) )
!!$                      imj = Latt%imj(i,j)
!!$                      Obs_eq(6)%Obs_Latt(imj,1,1,1) = Obs_eq(6)%Obs_Latt(imj,1,1,1) + &
!!$                           & cmplx(X,0.d0,kind(0.d0))* ZP * ZS
!!$                   Enddo
!!$                   if (ny /= L2/2) then
!!$                      J = Latt%nnlist(J,0,1)
!!$                      X   = X*DW_Ising_tau( nsigma(L_bond(J ,1),ntau)*nsigma(L_bond(J ,1),ntau1) )
!!$                      imj = Latt%imj(i,j)
!!$                      Obs_eq(6)%Obs_Latt(imj,1,1,1) = Obs_eq(6)%Obs_Latt(imj,1,1,1) + &
!!$                           & cmplx(X,0.d0,kind(0.d0))* ZP * ZS
!!$                   Endif
!!$                   ! J = I + 2*ax
!!$                Enddo
!!$             Enddo
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
                   Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + &
                        & (   GRC(I1,J1,1) * GR(I1,J1,1) +  GRC(I1,J1,2) * GR(I1,J1,2)    + &
                        &    (GRC(I1,I1,2) - GRC(I1,I1,1))*(GRC(J1,J1,2) - GRC(J1,J1,1))    ) * ZP*ZS
                   ! SpinXY
                   ! c^d_(i,u) c_(i,d) c^d_(j,d) c_(j,u)  +  c^d_(i,d) c_(i,u) c^d_(j,u) c_(j,d)
                   Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + &
                     & (   GRC(I1,J1,1) * GR(I1,J1,2) +  GRC(I1,J1,2) * GR(I1,J1,1)    ) * ZP*ZS
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
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS, ZN, G(4,4)
          Integer :: IMJ, I, J, I1, J1, no_I, no_J, NT1, NT1P1, no, no1

          NT1 = NT
          If (NT == 0 ) NT1 = LTROT
          NT1P1 = NT1 + 1
          if (NT1P1 > Ltrot ) NT1P1 = NT1P1 - Ltrot
          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          If (NT == 0 ) then 
             DO I = 1,Size(Obs_tau,1)
                Obs_tau(I)%N = Obs_tau(I)%N + 1
                Obs_tau(I)%Ave_sign = Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
             ENDDO
          endif
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising"  ) then 
             ZN =  cmplx(dble(N_SUN),0.d0, kind(0.D0))
             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   ! Green
                   Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        & +  ZN * GT0(I1,J1,1) * ZP* ZS
                   
                   ! SpinZ
                   Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        &      - ZN*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
                   
                   ! SpinXY
                   Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        &      - ZN*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
                   
                   ! Den
                   Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        & + ( ZN*ZN*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1))*       &
                        &         (cmplx(1.d0,0.d0,kind(0.d0)) - G00(J1,J1,1))  -     &
                        &     ZN * GT0(I1,J1,1)*G0T(J1,I1,1)                                ) * ZP * ZS
                Enddo
                Obs_tau(4)%Obs_Latt0(no_I) = Obs_tau(4)%Obs_Latt0(no_I) + &
                     &         ZN*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1)) * ZP * ZS
             Enddo
             Do I = 1,Latt%N
                do no = 1,N_coord
                   select  case(no)
                   case(1) 
                      I1 = Latt%nnlist(I,1,0)
                   case(2)
                      I1 = Latt%nnlist(I,0,1)
                   end select
                   Do J = 1,Latt%N
                      imj = latt%imj(I,J)
                      do no1 = 1,N_coord
                         if ( no1 == 1) J1 = Latt%nnlist(J,1,0)
                         if ( no1 == 2) J1 = Latt%nnlist(J,0,1)
                         G(1,1) = ZKRON(I,I ) - GTT(I ,I ,1)
                         G(1,2) = ZKRON(I,I1) - GTT(I1,I ,1)
                         G(1,3) = -G0T (J ,I,1)
                         G(1,4) = -G0T (J1,I,1)
                         G(2,1) = -GTT(I,I1 ,1)
                         G(2,2) =  ZKRON(I1,I1) -GTT(I1,I1,1)
                         G(2,3) = -G0T (J ,I1,1)
                         G(2,4) = -G0T (J1,I1,1)
                         G(3,1) = -GT0 (I ,J ,1)
                         G(3,2) = -GT0 (I1,J ,1)
                         G(3,3) =  ZKRON(J,J ) - G00(J ,J ,1)
                         G(3,4) =  ZKRON(J,J1) - G00(J1,J ,1)
                         G(4,1) = -GT0(I ,J1,1)
                         G(4,2) = -GT0(I1,J1,1)
                         G(4,3) = -G00(J ,J1,1)
                         G(4,4) = ZKRON(J1,J1) - G00(J1,J1,1)
!!$                         G(1,1)=  GRC(I,I,1); G(1,2)= GRC(I,I1,1) ;G(1,3)= GRC(I,J,1) ;G(1,4)=GRC(I,J1,1)
!!$                         G(2,1)= -GR(I,I1,1); G(2,2)= GRC(I1,I1,1);G(2,3)= GRC(I1,J,1);G(2,4)=GRC(I1,J1,1)
!!$                         G(3,1)= -GR(I,J,1) ; G(3,2)=-GR(I1,J,1)  ;G(3,3)= GRC(J,J,1) ;G(3,4)=GRC(J,J1,1)
!!$                         G(4,1)= -GR(I,J1,1); G(4,2)=-GR(I1,J1,1) ;G(4,3)=-GR(J,J1,1) ;G(4,4)=GRC(J1,J1,1)
                         If ( N_SUN == 4 ) then
#include "Bid_N4.f90"
                         endif
                         If ( N_SUN == 3 ) then
#include "Bid_N3.f90"
                         endif
                         If ( N_SUN == 2 ) then
#include "Bid_N2.f90"
                         endif
                         Obs_tau(7)%Obs_Latt(imj,nt+1,no,no1)  = Obs_tau(7)%Obs_Latt(imj,nt+1,no,no1)   +  Z * ZP*ZS

                         !Kin-Kin correlations
                         Z =  (  (GTT(I,I1,1) +  GTT(I1,I,1)) * (G00(J,J1,1)  +  G00(J1,J,1)) * ZN  - &
                              &   G0T(J1,I ,1)*GT0(I1,J ,1) - &
                              &   G0T(J ,I ,1)*GT0(I1,J1,1) - &
                              &   G0T(J1,I1,1)*GT0(I ,J ,1) - &
                              &   G0T(J ,I1,1)*GT0(I ,J1,1)   )* ZN * &
                              &   real(nsigma(L_bond(I ,no),NT1)*nsigma(L_bond(J ,no1),Ltrot), Kind(0.d0))
                         Obs_tau(6)%Obs_Latt(imj,nt+1,no,no1)  = Obs_tau(6)%Obs_Latt(imj,nt+1,no,no1)   +  Z * ZP*ZS

                         Z =  cmplx(DW_Ising_tau( nsigma(L_bond(J,no1),Ltrot)*nsigma(L_bond(J,no1),1    ) ) *&
                              &     DW_Ising_tau( nsigma(L_bond(I,no ),NT1  )*nsigma(L_bond(I,no ),NT1P1) ), 0.d0,kind(0.d0))
                         Obs_tau(5)%Obs_Latt(imj,nt+1,no,no1)  = Obs_tau(5)%Obs_Latt(imj,nt+1,no,no1)   +  Z * ZP*ZS
                      enddo
                      
                   enddo
                   If ( N_SUN == 2 ) then
#include "Bid_N2_0.f90"
                   endif
                   Obs_tau(7)%Obs_Latt0(no) =  Obs_tau(7)%Obs_Latt0(no) +  Z * ZP*ZS
                   Z = -(GTT(I,I1,1) +  GTT(I1,I,1)) *ZN * Real(nsigma(L_bond(I ,no),NT1),kind(0.d0))
                   Obs_tau(6)%Obs_Latt0(no) =  Obs_tau(6)%Obs_Latt0(no) +  Z * ZP*ZS
                   Z = cmplx( DW_Ising_tau( nsigma(L_bond(I,no ),NT1  )*nsigma(L_bond(I,no ),NT1P1) ), 0.d0,kind(0.d0))
                   Obs_tau(5)%Obs_Latt0(no)  = Obs_tau(5)%Obs_Latt0(no) +  Z * ZP*ZS
                enddo
             enddo
             
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

!===================================================================================           

      Integer Function  iFlux(I,nt)

        Implicit none
        
        Integer, INTENT(IN) :: I,nt

        ! Local 
        Integer :: n1,n2,n3,n4

        !   I3  I2
        !   I   I1
        n1  = L_bond(I,1)
        n2  = L_bond(Latt%nnlist(I,1,0),2)
        n3  = L_bond(Latt%nnlist(I,0,1),1)
        n4  = L_bond(I,2)
        iFlux =   nsigma(n1,nt)*nsigma(n2,nt)*nsigma(n3,nt)*nsigma(n4,nt)
        
      end Function iFlux

!===================================================================================
      Subroutine  Hamiltonian_set_random_nsigma

        ! The user can set the initial configuration
        
        Implicit none
        
        Integer :: I,nt, n

        Do nt = 1,Ltrot
           Do I = 1, Latt%N
              if (mod( Latt%list(i,1) + latt%list(i,2), 2 ) == 0 ) then 
                 nsigma(L_bond(I,1),nt) =  1
                 nsigma(L_bond(I,2),nt) = -1
              else
                 nsigma(L_bond(I,1),nt) =  1
                 nsigma(L_bond(I,2),nt) = 1
              endif
           Enddo
           Do n = N_coord*Latt%N +  1, Size(OP_V,1) ! This is needed in case U = 0
              nsigma(n,nt) = 1
              if ( ranf_wrap() > 0.5 ) nsigma(n,nt) = -1
           enddo
        enddo
        
        Call Print_fluxes
        
        
      end Subroutine Hamiltonian_set_random_nsigma
!===================================================================================

      Subroutine Print_fluxes

        Implicit none

#if defined (MPI) || defined(TEMPERING)
        include 'mpif.h'

#endif   
        

        ! Local 
        Integer :: I,nt,ix, iy, n
        Character (len=64) :: File1


#ifdef MPI
        Integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
        Integer        :: STATUS(MPI_STATUS_SIZE)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
        call MPI_Comm_rank(Group_Comm, irank_g, ierr)
        call MPI_Comm_size(Group_Comm, isize_g, ierr)
        igroup           = irank/isize_g
#endif

#if defined(TEMPERING)
        write(File1,'(A,I0,A)') "Temp_",igroup,"/Fluxes"
#else
        File1="Fluxes"
#endif

        Open (Unit=10,File=File1, status="unknown")
        Do nt = 1,Ltrot
           Do i  = 1,Ndim
              n = iFlux(I,nt)
              if (n == 1 ) then
                 ix = Latt%list(i,1)
                 iy = Latt%list(i,2)
                 Write(10,'(I4,2x,I4,2x,I4)')   IX, IY, NT
              endif
           Enddo
        enddo
        close(10)
      end Subroutine Print_fluxes
!===================================================================================


      Subroutine Test_Hamiltonian
        
        Implicit none
        
        Integer :: n,  nc, n_op, nt
        Integer, allocatable :: nsigma_old(:,:)
        Real (Kind=kind(0.d0)) :: X, X1, size_clust
        
        n = size(Op_V,1)
        allocate (nsigma_old(n,Ltrot))
        do nc = 1,100
           !nt  = nranf(Ltrot)
           !n_op= nranf(n)
           !if ( OP_V(n_op,1)%type == 1 ) then
           !   X = S0(n_op,nt)
           !   nsigma_old = nsigma
           !   nsigma(n_op,nt) = -nsigma(n_op,nt)
           !   X1 = Delta_S0_global(Nsigma_old) 
           !   Write(6,*) nc, X, X1
           !endif
           nsigma_old = nsigma
           Call Global_move(X,nsigma_old,size_clust)
           X1 = Delta_S0_global(Nsigma_old) 
           Write(6,*) nc, X, X1
        enddo
        deallocate (nsigma_old)

        stop

      end Subroutine Test_Hamiltonian

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


        end Subroutine Global_move_tau

      end Module Hamiltonian

