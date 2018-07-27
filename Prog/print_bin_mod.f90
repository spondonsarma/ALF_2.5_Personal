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

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Handles printing of bins etc. This module is outdated and replaces by the Observables module
!
!--------------------------------------------------------------------
    Module Print_bin_mod

       Interface Print_bin
          module procedure Print_bin_C, Print_bin_R
       end Interface Print_bin


       Contains

         Subroutine  Print_bin_C(Dat_eq,Dat_eq0,Latt, Nobs, Phase_bin_tmp, file_pr)
           Use Lattices_v3
#ifdef MPI
           Use mpi
#endif
           Implicit none


           Complex (Kind=Kind(0.d0)), Dimension(:,:,:), Intent(inout):: Dat_eq
           Complex (Kind=Kind(0.d0)), Dimension(:)    , Intent(inout):: Dat_eq0
           Type (Lattice),                     Intent(In)   :: Latt
           Complex (Kind=Kind(0.d0)),                   Intent(In)   :: Phase_bin_tmp
           Character (len=64),                 Intent(In)   :: File_pr
           Integer,                            Intent(In)   :: Nobs

           ! Local
           Integer :: Norb, I, no, no1
           Complex (Kind=Kind(0.d0)), allocatable :: Tmp(:,:,:)
           Real    (Kind=Kind(0.d0))              :: x_p(2) 
           Complex (Kind=Kind(0.d0))              :: Phase_bin
#ifdef MPI
           Complex (Kind=Kind(0.d0)), allocatable :: Tmp1(:)
           Complex (Kind=Kind(0.d0)):: Z
           Integer         :: Ierr, Isize, Irank
           INTEGER         :: STATUS(MPI_STATUS_SIZE)
           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
           Phase_bin = Phase_bin_tmp
           Norb = size(Dat_eq,3)
           if ( .not. (Latt%N  == Size(Dat_eq,1) ) ) then 
              Write(6,*) 'Error in Print_bin' 
              Stop
           endif
           Allocate (Tmp(Latt%N,Norb,Norb))
           Dat_eq = Dat_eq/dble(Nobs)
           Dat_eq0 = Dat_eq0/dble(Nobs*Latt%N)

#ifdef MPI
           I = Latt%N*Norb*Norb
           Tmp = cmplx(0.d0, 0.d0, kind(0.D0))
           CALL MPI_REDUCE(Dat_eq,Tmp,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Dat_eq = Tmp/DBLE(ISIZE)
           I = 1
           Z = cmplx(0.d0,0.d0,kind(0.d0))
           CALL MPI_REDUCE(Phase_bin,Z,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Phase_bin= Z/DBLE(ISIZE)

           I = Norb
           Allocate (Tmp1(Norb))
           Tmp1 = cmplx(0.d0,0.d0,kind(0.d0))
           CALL MPI_REDUCE(Dat_eq0,Tmp1,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Dat_eq0 = Tmp1/DBLE(ISIZE)
           deallocate(Tmp1)
           If (Irank == 0 ) then
#endif
              do no = 1,Norb
                 do no1 = 1,Norb
                    Call  Fourier_R_to_K(Dat_eq(:,no,no1), Tmp(:,no,no1), Latt)
                 enddo
              enddo
              Open (Unit=10,File=File_pr, status="unknown",  position="append")
              Write(10,*) dble(Phase_bin),Norb,Latt%N
              do no = 1,Norb
                 Write(10,*) Dat_eq0(no)
              enddo
              do I = 1,Latt%N
                 x_p = dble(Latt%listk(i,1))*Latt%b1_p + dble(Latt%listk(i,2))*Latt%b2_p  
                 Write(10,*) X_p(1), X_p(2)
                 do no = 1,Norb
                    do no1 = 1,Norb
                       Write(10,*) tmp(I,no,no1)
                    enddo
                 enddo
              enddo
              close(10)
#ifdef MPI
           Endif
#endif

           deallocate (Tmp)

         End Subroutine Print_bin_C

!=========================================================

         Subroutine  Print_bin_R(Dat_eq,Dat_eq0,Latt, Nobs, Phase_bin_tmp, file_pr)
           Use Lattices_v3
#ifdef MPI
           Use mpi
#endif
           Implicit none

           Real    (Kind=Kind(0.d0)), Dimension(:,:,:), Intent(inout) :: Dat_eq
           Real    (Kind=Kind(0.d0)), Dimension(:)    , Intent(inout) :: Dat_eq0
           Type (Lattice),                     Intent(In)    :: Latt
           Complex (Kind=Kind(0.d0)),                   Intent(In)    :: Phase_bin_tmp
           Character (len=64),                 Intent(In)    :: File_pr
           Integer,                            Intent(In)    :: Nobs
           
           ! Local
           Integer :: Norb, I, no,no1
           Real    (Kind=Kind(0.d0)), allocatable :: Tmp(:,:,:)
           Real    (Kind=Kind(0.d0))              :: x_p(2) 
           Complex (Kind=Kind(0.d0))              :: Phase_bin
#ifdef MPI
           Integer        :: Ierr, Isize, Irank
           Complex (Kind=Kind(0.d0))              :: Z
           Real    (Kind=Kind(0.d0)), allocatable :: Tmp1(:)
           INTEGER        :: STATUS(MPI_STATUS_SIZE)
           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
           
           Phase_bin = Phase_bin_tmp
           Norb = size(Dat_eq,3)
           if ( .not. (Latt%N  == Size(Dat_eq,1) ) ) then 
              Write(6,*) 'Error in Print_bin' 
              Stop
           endif
           Allocate (Tmp(Latt%N,Norb,Norb))
           Dat_eq  = Dat_eq/dble(Nobs)
           Dat_eq0 = Dat_eq0/(dble(Nobs)*dble(Latt%N))
#ifdef MPI
           I = Latt%N*Norb*Norb
           Tmp = 0.d0
           CALL MPI_REDUCE(Dat_eq,Tmp,I,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Dat_eq = Tmp/DBLE(ISIZE)
           I = 1
           Z = cmplx(0.d0,0.d0,kind(0.d0))
           CALL MPI_REDUCE(Phase_bin,Z,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Phase_bin= Z/DBLE(ISIZE)
           If (Irank == 0 ) then

           I = Norb
           Allocate (Tmp1(Norb) )
           Tmp1 = 0.D0
           CALL MPI_REDUCE(Dat_eq0,Tmp1,I,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Dat_eq0 = Tmp1/DBLE(ISIZE)
           deallocate( Tmp1)
#endif
              do no = 1,Norb
                 do no1 = 1,Norb
                    Call  Fourier_R_to_K(Dat_eq(:,no,no1), Tmp(:,no,no1), Latt)
                 enddo
              enddo
              Open (Unit=10,File=File_pr, status="unknown",  position="append")
              Write(10,*) dble(Phase_bin),Norb,Latt%N
              do no = 1,Norb
                 Write(10,*) Dat_eq0(no)
              enddo
              do I = 1,Latt%N
                 x_p = dble(Latt%listk(i,1))*Latt%b1_p + dble(Latt%listk(i,2))*Latt%b2_p  
                 Write(10,*) X_p(1), X_p(2)
                 do no = 1,Norb
                    do no1 = 1,Norb
                       Write(10,*) tmp(I,no,no1)
                    enddo
                 enddo
              enddo
              close(10)
#ifdef MPI
           endif
#endif
           deallocate (Tmp )
           
         End Subroutine Print_bin_R
!============================================================
         Subroutine  Print_scal(Obs, Nobs, file_pr)
#ifdef MPI
           Use mpi
#endif
           Implicit none

           Complex   (Kind=Kind(0.d0)), Dimension(:), Intent(inout) :: Obs
           Character (len=64),               Intent(In)    :: File_pr
           Integer,                          Intent(In)    :: Nobs

           ! Local
           Integer :: I
#ifdef MPI
           Complex  (Kind=Kind(0.d0)), allocatable :: Tmp(:)
           Integer        :: Ierr, Isize, Irank, Norb
           INTEGER        :: STATUS(MPI_STATUS_SIZE)
           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
           Obs = Obs/dble(Nobs)
#ifdef MPI
           Norb = size(Obs,1)
           Allocate ( Tmp(Norb) )
           Tmp = 0.d0
           CALL MPI_REDUCE(Obs,Tmp,Norb,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Obs = Tmp/DBLE(ISIZE)
           if (Irank == 0 ) then
#endif
              Open (Unit=10,File=File_pr, status="unknown",  position="append")
              WRITE(10,*) (Obs(I), I=1,size(Obs,1))
              close(10)
#ifdef MPI
           endif
           deallocate (Tmp )
#endif
           
         End Subroutine Print_scal

!==============================================================
         Subroutine  Print_bin_tau(Dat_tau, Latt, Nobs, Phase_bin, file_pr, dtau, Dat0_tau)
           Use Lattices_v3
#ifdef MPI
           Use mpi
#endif
           Implicit none

           Complex (Kind=Kind(0.d0)), Dimension(:,:,:,:), Intent(inout):: Dat_tau   ! (Latt%N, Ltau,Norb, Norb)
           Complex (Kind=Kind(0.d0)), Dimension(:      ), Intent(inout), optional :: Dat0_tau  ! (Norb)
           Type (Lattice),                       Intent(In)   :: Latt
           Complex (Kind=Kind(0.d0)),                     Intent(In)   :: Phase_bin
           Character (len=64),                   Intent(In)   :: File_pr
           Integer,                              Intent(In)   :: Nobs
           Real (Kind=Kind(0.d0)),                        Intent(In)   :: dtau

           ! Local
           Integer :: Norb, I, no,no1, LT, nt
           Complex (Kind=Kind(0.d0)), allocatable :: Tmp(:,:,:,:), Tmp0(:)
           Complex (Kind=Kind(0.d0)) :: Phase_mean 
           Real    (Kind=Kind(0.d0))              :: x_p(2) 
#ifdef MPI
           Complex (Kind=Kind(0.d0)):: Z
           Integer         :: Ierr, Isize, Irank
           INTEGER         :: STATUS(MPI_STATUS_SIZE)
           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif

           Phase_mean = Phase_bin
           Norb = size(Dat_tau,3)
           if ( .not. (Latt%N  == Size(Dat_tau,1) ) ) then 
              Write(6,*) 'Error in Print_bin' 
              Stop
           endif
           LT = Size(Dat_tau,2)
           Allocate (Tmp (Latt%N,LT,Norb,Norb) )
           Allocate (Tmp0(Norb) )

           Tmp0 = cmplx(0.d0, 0.d0, kind(0.D0))
           Dat_tau  = Dat_tau/dble(Nobs)
           If (Present(Dat0_tau) ) Dat0_tau = Dat0_tau/dble(Nobs*Latt%N*LT)
           
#ifdef MPI
           I = Latt%N*Norb*Norb*LT
           Tmp = cmplx(0.d0, 0.d0, kind(0.D0))
           CALL MPI_REDUCE(Dat_tau,Tmp,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Dat_tau = Tmp/DBLE(ISIZE)
           
           If (Present(Dat0_tau) ) then
              I = Norb
              Tmp0 = cmplx(0.d0, 0.d0, kind(0.D0))
              CALL MPI_REDUCE(Dat0_tau,Tmp0,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
              Dat0_tau = Tmp0/DBLE(ISIZE)
           endif

           I = 1
           Z = cmplx(0.d0,0.d0,kind(0.d0))
           CALL MPI_REDUCE(Phase_mean,Z,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Phase_mean= Z/DBLE(ISIZE)
           If (Irank == 0 ) then
#endif
              If (Present(Dat0_tau) ) Tmp0 = Dat0_tau
              do nt = 1,LT
                 do no = 1,Norb
                    do no1 = 1,Norb
                       Call  Fourier_R_to_K(Dat_tau(:,nt,no,no1), Tmp(:,nt,no,no1), Latt)
                    enddo
                 enddo
              enddo
              Open (Unit=10,File=File_pr, status="unknown",  position="append")
              Write(10,*) dble(Phase_mean),Norb,Latt%N, LT, dtau
              Do no = 1,Norb
                 Write(10,*)  Tmp0(no)
              enddo
              do I = 1,Latt%N
                 x_p = dble(Latt%listk(i,1))*Latt%b1_p + dble(Latt%listk(i,2))*Latt%b2_p  
                 Write(10,*) X_p(1), X_p(2)
                 Do nt = 1,LT
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
              
           deallocate (Tmp, Tmp0 )

         End Subroutine Print_bin_tau


         
       end Module Print_bin_mod
