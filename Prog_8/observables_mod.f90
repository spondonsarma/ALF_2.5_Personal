     Module Observables

       
!--------------------------------------------------------------------
!> @author
!> Fakher Assaad
!
!> @brief 
!> This module. 
!> 1) Obser_Vec 
!> 2) Obser_Latt
!--------------------------------------------------------------------
    
       Type Obser_Vec
!>  Data structure for 
!>  < O_n >  n : =1, size(Obs,1)  
          Integer            :: N                    ! Number of measurements
          complex   (kind=8) :: Phase                ! Averarge phase
          complex   (kind=8), pointer :: Obs_vec(:)  ! Vector of observables
          Character (len=64) :: File_Vec             ! Name of file in which the bins will be written out
       end type Obser_Vec
       

       Type Obser_Latt
!>  Data structure for 
!>  < O^{dagger}(i,tau)_n O(j,0)_m>  - < O_n> <O_m> 
!>  where it is assumed that translation symmetry as specified by the lattice Latt is present.  
!>  Obs_Latt(i-j,tau,n,m) = < O^{dagger}(i,tau)_n O(j,0)_m>  
!>  Obs_Latt0(n) = < O_n>
!>  For equal   time correlation functions, tau runs from 1,1 
!>  For unequal time correlation functions, tau runs from 1,Ltrot+1  
          Integer            :: N                           ! Number of measurements
          complex   (kind=8) :: Phase                       ! Averarge phase
          complex   (kind=8), pointer :: Obs_Latt (:,:,:,:) ! i-j, tau, norb, norb  
          complex   (kind=8), pointer :: Obs_Latt0(:)       ! norb 
          Character (len=64) :: File_Latt                   ! Name of file in which the bins will be written out
       end type Obser_Latt
       


       Contains

         Subroutine Obser_Latt_make(Obs,Ns,Nt,No,Filename)
           Implicit none
           Type (Obser_Latt), intent(INOUT) :: Obs
           Integer, Intent(IN)             :: Ns,Nt,No
           Character (len=64), Intent(IN)  :: Filename 
           Allocate (Obs%Obs_Latt (Ns,Nt,No,No))
           Allocate (Obs%Obs_Latt0(No)         )
           Obs%File_Latt = Filename
         end subroutine Obser_Latt_make

         Subroutine Obser_Latt_Init(Obs)
           Implicit none
           Type (Obser_Latt), intent(INOUT) :: Obs
           Obs%Obs_Latt  = cmplx(0.d0,0.d0,kind(0.d0))
           Obs%Obs_Latt0 = cmplx(0.d0,0.d0,kind(0.d0))
           Obs%N         = 0
           Obs%Phase     = cmplx(0.d0,0.d0,kind(0.d0))
         end subroutine Obser_Latt_Init
         

         Subroutine Obser_Vec_make(Obs,N,Filename)
           Implicit none
           Type (Obser_vec), intent(INOUT) :: Obs
           Integer, Intent(IN)             :: N
           Character (len=64), Intent(IN)  :: Filename 
           Allocate (Obs%Obs_vec(N))
           Obs%File_Vec = Filename
         end subroutine Obser_Vec_make
         
         Subroutine Obser_Vec_Init(Obs)
           Implicit none
           Type (Obser_vec), intent(INOUT) :: Obs
           Obs%Obs_vec = cmplx(0.d0,0.d0,kind(0.d0))
           Obs%N       = 0
           Obs%Phase   = cmplx(0.d0,0.d0,kind(0.d0))
         end subroutine Obser_Vec_Init

         

!!!!!!!!!!!!!!!         
         
         Subroutine  Print_bin_Latt(Obs,Latt,dtau)
           Use Lattices_v3
           Implicit none

#include "machine"
#ifdef MPI
           include 'mpif.h'
#endif   

           Type (Obser_Latt),        Intent(Inout)   :: Obs
           Type (Lattice),           Intent(In)      :: Latt
           Real (Kind=8),            Intent(In)      :: dtau

           ! Local
           Integer :: Ns,Nt, Norb, no, no1, I 
           Complex (Kind=8), allocatable :: Tmp(:,:,:,:), Tmp1(:)
           Real    (Kind=8)              :: x_p(2) 
           Complex (Kind=8)              :: Phase_bin

#ifdef MPI
           Complex (Kind=8):: Z
           Integer         :: Ierr, Isize, Irank
           INTEGER         :: STATUS(MPI_STATUS_SIZE)
           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif

           Ns    = Size(Obs%Obs_Latt,1)
           Nt    = Size(Obs%Obs_Latt,2)
           Norb  = Size(Obs%Obs_Latt,3)
           if ( .not. (Latt%N  == Ns ) ) then 
              Write(6,*) 'Error in Print_bin' 
              Stop
           endif
           Allocate (Tmp(Ns,Nt,Norb,Norb), Tmp1(Norb) )
           Obs%Obs_Latt  =   Obs%Obs_Latt /dble(Obs%N   )
           Obs%Obs_Latt0 =   Obs%Obs_Latt0/dble(Obs%N*Ns)
           Obs%Phase     =   Obs%Phase    /dble(Obs%N   )

#ifdef MPI
           I = Ns*Nt*Norb*Norb
           Tmp = cmplx(0.d0, 0.d0, kind(0.D0))
           CALL MPI_REDUCE(Obs%Obs_Latt,Tmp,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Obs%Obs_Latt = Tmp/DBLE(ISIZE)

           I = 1
           Z = cmplx(0.d0,0.d0,kind(0.d0))
           CALL MPI_REDUCE(Obs%Phase,Z,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Obs%Phase = Z/DBLE(ISIZE)

           I = Norb
           Tmp1 = cmplx(0.d0,0.d0,kind(0.d0))
           CALL MPI_REDUCE(Obs%Obs_Latt0,Tmp1,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Obs%Obs_Latt0 = Tmp1/DBLE(ISIZE)

           If (Irank == 0 ) then
#endif
              do nt = 1,Nt
                 do no = 1,Norb
                    do no1 = 1,Norb
                       Call  Fourier_R_to_K(Obs%Obs_Latt(:,nt,no,no1), Tmp(:,nt,no,no1), Latt)
                    enddo
                 enddo
              enddo
              Open (Unit=10,File=Obs%File_Latt, status="unknown",  position="append")
              Write(10,*) dble(Obs%Phase),Norb,Latt%N, Nt, dtau
              Do no = 1,Norb
                 Write(10,*)  Obs%Obs_Latt0(no)
              enddo
              do I = 1,Latt%N
                 x_p = dble(Latt%listk(i,1))*Latt%b1_p + dble(Latt%listk(i,2))*Latt%b2_p  
                 Write(10,*) X_p(1), X_p(2)
                 Do nt = 1,Nt
                    do no = 1,Norb
                       do no1 = 1,Norb
                          Write(10,*) tmp(I,nt,no,no1)
                       enddo
                    enddo
                 enddo
              enddo
              close(10)
#ifdef MPI
           Endif
#endif
              
           deallocate (Tmp, tmp1 )
          

         End Subroutine Print_bin_Latt
         

!============================================================
         Subroutine  Print_bin_Vec(Obs)
           
           Implicit none

#include "machine"
#ifdef MPI
           include 'mpif.h'
#endif   
           
           Type (Obser_vec), intent(Inout) :: Obs
           
           ! Local
           Integer :: No,I
           Complex  (Kind=8), allocatable :: Tmp(:)
#ifdef MPI
           Integer        :: Ierr, Isize, Irank
           INTEGER        :: STATUS(MPI_STATUS_SIZE)
           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
           
           No = size(Obs%Obs_vec,1)
           Allocate (Tmp(No) )
           Obs%Obs_vec = Obs%Obs_vec/dble(Obs%N)
           Obs%Phase   = Obs%Phase  /dble(Obs%N)


#ifdef MPI
           Tmp = cmplx(0.d0,0.d0,kind(0.d0))
           CALL MPI_REDUCE(Obs%Obs_vec,Tmp,No,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Obs%Obs_vec = Tmp/DBLE(ISIZE)
           if (Irank == 0 ) then
#endif
              Open (Unit=10,File=Obs%File_Vec, status="unknown",  position="append")
              WRITE(10,*) (Obs%Obs_vec(I), I=1,size(Obs%Obs_vec,1)), Obs%Phase
              close(10)
#ifdef MPI
           endif
#endif
           deallocate (Tmp )
           
         End Subroutine Print_bin_Vec

         
       end Module Observables
