    !Model Hamiltonian for interaction-induced topological reduction
    Module Hamiltonian

      Use Operator_mod
      Use Lattices_v3 
      Use MyMats 
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      

      Type (Operator), dimension(:,:), allocatable  :: Op_V
      Type (Operator), dimension(:,:), allocatable  :: Op_T
      Integer, allocatable :: nsigma(:,:)
      Integer              :: Ndim,  N_FL,  N_SUN,  Ltrot


      
      ! What is below is  private 
      
      Type (Lattice),       private :: Latt
      Integer, parameter,   private :: Norb=16
      Integer, allocatable, private :: List(:,:), Invlist(:,:)
      Integer,              private :: L1, L2
      real (Kind=8),        private :: Ham_T, Ham_Vint,  Ham_Lam
      real (Kind=8),        private :: Dtau, Beta
      Character (len=64),   private :: Model, Lattice_type
      Complex (Kind=8),     private :: Gamma_M(4,4,5), Sigma_M(2,2,0:3)


      ! Observables
      Integer,                       private :: Nobs
      Complex (Kind=8), allocatable, private :: obs_scal(:)
      Complex (Kind=8), allocatable, private ::  Den_eq(:,:,:), Den_eq0(:)

      ! For time displaced
      Integer,                       private :: NobsT
      Complex (Kind=8),              private :: Phase_tau
      Complex (Kind=8), allocatable, private :: Green_tau(:,:,:,:), Den_tau(:,:,:,:)

      contains 

        Subroutine Ham_Set

          Implicit none
#include "machine"
#ifdef MPI
          include 'mpif.h'
#endif   

          integer :: ierr

          NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model

          NAMELIST /VAR_SPT/  ham_T, Ham_Vint,  Ham_Lam,  Dtau, Beta


#ifdef MPI
          Integer        :: Isize, Irank
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
          
          
          OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
          IF (ierr /= 0) THEN
             WRITE(*,*) 'unable to open <parameters>',ierr
             STOP
          END IF
          READ(5,NML=VAR_lattice)
          CLOSE(5)
#ifdef MPI
          CALL MPI_BCAST(L1          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(L2          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
#endif
          Call Ham_latt

          N_FL  = 1
          N_SUN = 1
          
#ifdef MPI
          If (Irank == 0 ) then
#endif
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             READ(5,NML=VAR_SPT)
             CLOSE(5)
#ifdef MPI
          endif

          CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ham_V    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ham_Lam  ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Dtau     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Beta     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif

          Call Ham_hop
          Ltrot = nint(beta/dtau)
#ifdef MPI
          If (Irank == 0) then
#endif
             Open (Unit = 50,file="info",status="unknown",position="append")
             Write(50,*) '====================================='
             Write(50,*) 'Model is      : ', Model 
             Write(50,*) 'Lattice is    : ', Lattice_type
             Write(50,*) 'Beta          : ', Beta
             Write(50,*) 'dtau,Ltrot    : ', dtau,Ltrot
             Write(50,*) 't             : ', Ham_T
             Write(50,*) 'V             : ', Ham_Vint
             Write(50,*) 'Lambda        : ', Ham_Lam
             close(50)
#ifdef MPI
          endif
#endif
          call Ham_V
        end Subroutine Ham_Set
!=============================================================================

        Subroutine Ham_Latt
          Implicit none
          !Set the lattice
          Integer :: no, I, nc
          Real (Kind=8)  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          If ( Lattice_type =="Square" ) then
             a1_p(1) =  1.0  ; a1_p(2) =  0.d0
             a2_p(1) =  0.0  ; a2_p(2) =  1.d0
             L1_p    =  dble(L1)*a1_p
             L2_p    =  dble(L2)*a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
             !Write(6,*)  'Lattice: ', Ndim
          else
             Write(6,*) "Lattice not yet implemented!"
             Stop
          endif

          Ndim = Latt%N*Norb
          Allocate (List(Ndim,Norb), Invlist(Latt%N,Norb))
          nc = 0
          Do I = 1,Latt%N
             Do no = 1,Norb
                nc = nc + 1
                List(nc,1) = I
                List(nc,2) = no
                Invlist(I,no) = nc 
                ! no = 1..4   Xi_1
                ! no = 5..8   Xi_2
                ! no = 9..12  Xi_3
                ! no = 13..16 Xi_4
             Enddo
          Enddo

        end Subroutine Ham_Latt

!===================================================================================           
        Subroutine Ham_hop
          Implicit none

          !  Setup the hopping
          !  Per flavor, the  hopping is given by: 
          !  e^{-dtau H_t}  =    Prod_{n=1}^{Ncheck} e^{-dtau_n H_{n,t}}

          Integer :: I, I1, I2,I3,no, no1,  n, Ncheck, nc , nth
          Real    (Kind=8) :: X
          Complex (Kind=8) :: Z


          ! Setup Gamma matrices
          Gamma_M = cmplx(0.d0,0.d0)
          Sigma_M = cmplx(0.d0,0.d0)
          Sigma_M(1,1,0) = cmplx( 1.d0, 0.d0)
          Sigma_M(2,2,0) = cmplx( 1.d0, 0.d0)
          Sigma_M(1,2,1) = cmplx( 1.d0, 0.d0)
          Sigma_M(2,1,1) = cmplx( 1.d0, 0.d0)
          Sigma_M(1,2,2) = cmplx( 0.d0,-1.d0)
          Sigma_M(2,1,2) = cmplx( 0.d0, 1.d0)
          Sigma_M(1,1,3) = cmplx( 1.d0, 0.d0)
          Sigma_M(2,2,3) = cmplx(-1.d0, 0.d0)
          Do no = 1,2
             Do no1 = 1,2
                Gamma_M(no+2,no1  ,1)  =  Sigma_M(no,no1,0) 
                Gamma_M(no  ,no1+2,1)  =  Sigma_M(no,no1,0) 
                Gamma_M(no+2,no1  ,2)  =  cmplx( 0.d0,-1.d0)*Sigma_M(no,no1,0)
                Gamma_M(no  ,no1+2,2)  =  cmplx( 0.d0, 1.d0)*Sigma_M(no,no1,0)
                Gamma_M(no  ,no1  ,3)  =  cmplx( 1.d0, 0.d0)*Sigma_M(no,no1,1)
                Gamma_M(no+2,no1+2,3)  =  cmplx(-1.d0, 0.d0)*Sigma_M(no,no1,1)
                Gamma_M(no  ,no1  ,4)  =  cmplx( 1.d0, 0.d0)*Sigma_M(no,no1,2)
                Gamma_M(no+2,no1+2,4)  =  cmplx(-1.d0, 0.d0)*Sigma_M(no,no1,2)
                Gamma_M(no  ,no1  ,5)  =  cmplx( 1.d0, 0.d0)*Sigma_M(no,no1,3)
                Gamma_M(no+2,no1+2,5)  =  cmplx(-1.d0, 0.d0)*Sigma_M(no,no1,3)
             Enddo
          Enddo

          Ncheck = 1
          allocate(Op_T(Ncheck,N_FL))
          do n = 1,N_FL
             Do nc = 1,NCheck
                Call Op_make(Op_T(nc,n),Ndim)
                DO I = 1, Ndim
                   Op_T(nc,n)%P(I) = I 
                enddo
                Do I = 1,Latt%N
                   do nth = 0,3
                      do no = 1,4
                         do no1 = 1,4
                            Z =  cmplx(1.d0 + 2.d0*Ham_Lam,0.d0)*Gamma_M(no,no1,3)
                            Op_T(nc,n)%O( invlist(I ,no  + 4*nth), invlist(I ,no1 + 4*nth ) )  =  Z
                         enddo
                      enddo
                      I1 =  Latt%nnlist(I,1,0)
                      do no = 1,4
                         do no1 = 1,4
                            Z =  (cmplx(0.d0,1.d0)*Gamma_M(no,no1,1) + Gamma_M(no,no1,3))/cmplx(2.d0,0.d0)
                            Op_T(nc,n)%O( invlist(I ,no  + 4*nth), invlist(I1,no1 + 4*nth ) )  =  Z
                            Op_T(nc,n)%O( invlist(I1,no1 + 4*nth), invlist(I ,no  + 4*nth ) )  =  conjg(Z)
                         enddo
                      enddo
                      I2   = Latt%nnlist(I,0,1)
                      do no = 1,4
                         do no1 = 1,4
                            Z =  (cmplx(0.d0,1.d0)*Gamma_M(no,no1,2) + Gamma_M(no,no1,3))/cmplx(2.d0,0.d0)
                            Op_T(nc,n)%O( invlist(I ,no  + 4*nth), invlist(I2,no1 + 4*nth ) )  =  Z
                            Op_T(nc,n)%O( invlist(I2,no1 + 4*nth), invlist(I ,no  + 4*nth ) )  =  conjg(Z)
                         enddo
                      enddo
                   enddo
                enddo
                Op_T(nc,n)%g=cmplx(-Dtau,0.d0)
                Call Op_set(Op_T(nc,n)) 
                ! Just for tests
                Do I = 1, Ndim
                   Write(6,*) Op_T(nc,n)%E(i)
                enddo
             enddo
          enddo


        end Subroutine Ham_hop
!===================================================================================           
        Subroutine Ham_V
          
          Implicit none 
          
          Integer :: nf, nth, n, n1, n2, n3, n4, I, I1, I2, J,  Ix, Iy, nc, no,no1, ns, npm 
          Integer :: nxy
          Real (Kind=8) :: X_p(2), X1_p(2), X2_p(2), X, XJ, Xpm

          Complex (Kind=8) :: Ps(4,4,2), Ps_G5(4,4,2), Tmp(4,4), Z
          Complex (Kind=8) :: Sx(16,16,2,2), Sy(16,16,2,2)


          Ps = cmplx(0.d0,0.d0)
          Call mmult (Tmp,Gamma_M(:,:,3), Gamma_M(:,:,4) )
          do ns = 1,2
             if (ns == 1) X =  1.d0/2d0
             if (ns == 2) X = -1.d0/2.d0
             Do I = 1,4
                Do J = 1,4
                   Z = cmplx(0.d0,0.d0)
                   if ( I == J )  Z = cmplx(1.d0/2.d0,0.d0)
                   Ps(I,J,ns) =   Z  + cmplx(0.d0,X) * tmp(I,J)
                Enddo
             Enddo
          Enddo
          
          Do ns = 1,2
             Call mmult ( Ps_G5(:,:,ns), Ps(:,:,ns), Gamma_M(:,:,5) )
          enddo
      
          Sx = cmplx(0.d0,0.d0)
          Sy = cmplx(0.d0,0.d0)
          Do ns = 1,2
             Do npm = 1,2
                if (npm == 1) Xpm =  1.0
                if (npm == 2) Xpm = -1.0
                Do no = 1,4
                   do no1 = 1,4
                      Sx(no    , no1 + 4 ,ns,npm) =  cmplx(1.d0, 0.d0)*Ps_G5(no,no1,ns)
                      Sx(no +4 , no1     ,ns,npm) =  cmplx(1.d0, 0.d0)*Ps_G5(no,no1,ns)
                      Sx(no +8 , no1 + 12,ns,npm) =  cmplx(xpm,  0.d0)*Ps_G5(no,no1,ns)
                      Sx(no+12 , no1 + 8 ,ns,npm) =  cmplx(xpm,  0.d0)*Ps_G5(no,no1,ns)
                      
                      Sy(no    , no1 + 4 ,ns,npm) =  cmplx(0.d0, -1.d0    )*Ps_G5(no,no1,ns)
                      Sy(no +4 , no1     ,ns,npm) =  cmplx(0.d0,  1.d0    )*Ps_G5(no,no1,ns)
                      Sy(no +8 , no1 + 12,ns,npm) =  cmplx(0.d0,  1.d0*xpm)*Ps_G5(no,no1,ns)
                      Sy(no+12 , no1 + 8 ,ns,npm) =  cmplx(0.d0, -1.d0*xpm)*Ps_G5(no,no1,ns)
                   enddo
                enddo
             enddo
          enddo


          ! Number of opertors 8 per unit cell
          Allocate( Op_V(8*Latt%N,N_FL) )
          do nf = 1,N_FL
             do i  = 1, 8*Latt%N
                Call Op_make(Op_V(i,nf),Norb) 
             enddo
          enddo
          nc = 0
          Do nf = 1,N_FL
             do nxy = 1,2
                do ns = 1,2
                   do npm = 1,2 
                      Xpm = 1.d0
                      if (npm == 2) Xpm = -1.d0
                      Do i = 1,Latt%N
                         nc = nc + 1 
                         Do no = 1,Norb
                            Op_V(nc,nf)%P(no)   = Invlist(I,no)  
                         enddo
                         Do no = 1,Norb
                            Do no1 = 1,Norb
                               If (nxy == 1)  Op_V(nc,nf)%O(no,no1) = Sx(no,no1,ns,npm)
                               If (nxy == 2)  Op_V(nc,nf)%O(no,no1) = Sy(no,no1,ns,npm)
                            Enddo
                         Enddo
                         Op_V(nc,nf)%g = SQRT(CMPLX(-Xpm*DTAU*Ham_Vint/8.d0,0.D0)) 
                         Op_V(nc,nf)%alpha  = cmplx(0.d0,0.d0)
                         Op_V(nc,nf)%type   = 2
                         Call Op_set( Op_V(nc,nf) )
                         ! The operator reads: 
                         ! g*s*( c^{dagger} O c   - alpha ))
                         ! with s the HS field.
                      Enddo
                   Enddo
                Enddo
             Enddo
          Enddo
          
        end Subroutine Ham_V
        
!===================================================================================           
        Real (Kind=8) function S0(n,nt)  
          Implicit none
          Integer, Intent(IN) :: n,nt 
          Integer :: i, nt1 
          S0 = 1.d0
        end function S0

!===================================================================================           
        Subroutine  Alloc_obs(Ltau) 

          Implicit none
          Integer, Intent(In) :: Ltau
          Integer :: I
          Allocate ( Obs_scal(5) )
          Allocate ( Den_eq(Latt%N,Norb,Norb), Den_eq0(Norb) ) 
          If (Ltau == 1) then 
             Allocate ( Green_tau(Latt%N,Ltrot+1,Norb,Norb), Den_tau(Latt%N,Ltrot+1,Norb,Norb) )
          endif
          
        end Subroutine Alloc_obs

!===================================================================================           
        
        Subroutine  Init_obs(Ltau) 

          Implicit none
          Integer, Intent(In) :: Ltau
          
          Integer :: I,n
          
          Nobs = 0
          Obs_scal  = cmplx(0.d0,0.d0)
          Den_eq    = cmplx(0.d0,0.d0)
          Den_eq0   = cmplx(0.d0,0.d0)

          If (Ltau == 1) then
             NobsT = 0
             Phase_tau = cmplx(0.d0,0.d0)
             Green_tau = cmplx(0.d0,0.d0)
             Den_tau = cmplx(0.d0,0.d0)
          endif

        end Subroutine Init_obs
        
!========================================================================
        Subroutine Obser(GR,Phase,Ntau)
          
          Implicit none
          
          Complex (Kind=8), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=8), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          
          !Local 
          Complex (Kind=8) :: GRC(Ndim,Ndim,N_FL), ZK
          Complex (Kind=8) :: Zrho, Zkin, ZPot, Z, ZP,ZS
          Integer :: I,J, no,no1, n, n1, imj, nf, I1, I2, J1, J2
          
          Real (Kind=8) :: G(4,4), X, FI, FJ
          
          Nobs = Nobs + 1
          ZP = PHASE/cmplx(Real(Phase,kind=8),0.d0)
          ZS = cmplx(Real(Phase,kind=8)/Abs(Real(Phase,kind=8)), 0.d0)
          

          Do nf = 1,N_FL
             Do I = 1,Ndim
                Do J = 1,Ndim
                   ZK = cmplx(0.d0,0.d0)
                   If ( I == J ) ZK = cmplx(1.d0,0.d0)
                   GRC(I,J,nf)  = ZK - GR(J,I,nf)
                Enddo
             Enddo
          Enddo
          ! GRC(i,j,nf) = < c^{dagger}_{j,nf } c_{j,nf } >
          ! Compute scalar observables. 

          Zkin = cmplx(0.d0,0.d0)
          Do nf = 1,N_FL
             Do J = 1,Op_T(1,nf)%N
                J1 = Op_T(1,nf)%P(J)
                DO I = 1,Op_T(1,nf)%N
                   I1 = Op_T(1,nf)%P(I)
                   Zkin  = Zkin  + Op_T(1,nf)%O(i,j)*Grc(i1,j1,nf) 
                Enddo
             ENddo
          Enddo
          Zkin = Zkin*cmplx( dble(N_SUN), 0.d0 )

          Zrho = cmplx(0.d0,0.d0)
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Zrho = Zrho + Grc(i,i,nf) 
             enddo
          enddo
          Zrho = Zrho*cmplx( dble(N_SUN), 0.d0 )

          ZPot = cmplx(0.d0,0.d0)

          Obs_scal(1) = Obs_scal(1) + zrho * ZP*ZS
          Obs_scal(2) = Obs_scal(2) + zkin * ZP*ZS
          Obs_scal(3) = Obs_scal(3) + Zpot * ZP*ZS
          Obs_scal(4) = Obs_scal(4) + (zkin +  Zpot)*ZP*ZS
          Obs_scal(5) = Obs_scal(5) + ZS
          ! You will have to allocate more space if you want to include more  scalar observables.
          DO I1 = 1,Ndim
             I  = List(I1,1)
             no = List(I1,2)
             DO J1 = 1, Ndim
                J = List(J1,1)
                no1 = list(J1,2)
                imj = latt%imj(I,J)

                DEN_Eq (imj,no,no1) = DEN_Eq (imj,no,no1)   +  &
                     & (   GRC(I1,J1,1) * GR (I1,J1,1)      +  &
                     &     GRC(I1,I1,1) * GRC(J1,J1,1)         ) * ZP*ZS

             enddo
             Den_eq0(no) = Den_eq0(no) +   GRC(I1,I1,1)*ZP*ZS 
          enddo

        end Subroutine Obser
!==========================================================        

        Subroutine  Pr_obs(LTAU)

          Use Print_bin_mod
          Implicit none
#include "machine"
#ifdef MPI
          include 'mpif.h'
#endif   


          Integer,  Intent(In) ::  Ltau

          Character (len=64) :: File_pr
          Complex   (Kind=8) :: Phase_bin
#ifdef MPI
          Integer        :: Isize, Irank, Ierr
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif          
!!$#ifdef MPI
!!$          Write(6,*)  Irank, 'In Pr_obs', LTAU
!!$#else
!!$          Write(6,*)  'In Pr_obs', LTAU
!!$#endif
    
          Phase_bin = Obs_scal(5)/cmplx(dble(Nobs),0.d0)
          File_pr ="Den_eq"
          Call Print_bin(Den_eq, Den_eq0, Latt, Nobs, Phase_bin, file_pr)

          File_pr ="ener"
          Call Print_scal(Obs_scal, Nobs, file_pr)
          If (Ltau == 1) then
             Phase_tau = Phase_tau/cmplx(dble(NobsT),0.d0)
             File_pr = "Green_tau"
             Call Print_bin_tau(Green_tau,Latt,NobsT,Phase_tau, file_pr,dtau)
             File_pr = "Den_tau"
             Call Print_bin_tau(Den_tau,Latt,NobsT,Phase_tau, file_pr,dtau)
          endif
!!$#ifdef MPI
!!$          Write(6,*)  Irank, 'out Pr_obs', LTAU
!!$#else
!!$          Write(6,*)  'out Pr_obs', LTAU
!!$#endif
        end Subroutine Pr_obs
!==========================================================        

        Subroutine OBSERT(NT,  GT0,G0T,G00,GTT, PHASE)
          Implicit none
          
          Integer         , INTENT(IN) :: NT
          Complex (Kind=8), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=8), INTENT(IN) :: Phase
          
          !Locals
          Complex (Kind=8) :: Z, ZP, ZS
          Integer :: IMJ, I, J

          ZP = PHASE/cmplx(Real(Phase,kind=8),0.d0)
          ZS = cmplx(Real(Phase,kind=8)/Abs(Real(Phase,kind=8)), 0.d0)
          If (NT == 0 ) then 
             Phase_tau = Phase_tau + ZS
             NobsT     = NobsT + 1
          endif
          If ( N_FL == 1 ) then 
             Z =  cmplx(dble(N_SUN),0.d0)
             Do I = 1,Latt%N
                Do J = 1,Latt%N
                   imj = latt%imj(I,J)
                   Green_tau(imj,nt+1,1,1) = green_tau(imj,nt+1,1,1)  +  Z * GT0(I,J,1) * ZP* ZS
                   Den_tau  (imj,nt+1,1,1) = Den_tau  (imj,nt+1,1,1)  -  Z * GT0(I,J,1)*G0T(J,I,1) * ZP* ZS
                Enddo
             Enddo
          Endif
        end Subroutine OBSERT


    end Module Hamiltonian
