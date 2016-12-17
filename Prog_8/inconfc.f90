!  Copyright (C) 2016 The ALF project
! 
!  This file is part of the ALF project.
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

        SUBROUTINE confin 

!--------------------------------------------------------------------
!
!> @brief 
!> Reads in the field configuration if present. If not present  the 
!> routine will generate one based on the seeds read in from the file
!> seeds
!
!--------------------------------------------------------------------
         Use Hamiltonian
         Implicit none
   
#ifdef MPI
         INCLUDE 'mpif.h'
         ! Local
#endif   

         Integer        :: I, IERR, ISIZE, IRANK, seed_in, K, iseed, Nt
         Integer, dimension(:), allocatable :: Seed_vec
         Real (Kind=8)  :: X
         Logical ::   lconf 
         character (len=64) :: file_sr, File_tg

#ifdef MPI
         INTEGER        :: STATUS(MPI_STATUS_SIZE)
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif

         Allocate (Nsigma(Size(Op_V,1),Ltrot))
         
#ifdef MPI
         INQUIRE (FILE='confin_0', EXIST=lconf)
         If (lconf) Then
            file_sr = "confin"
            Call Get_seed_Len(K)
            Allocate(Seed_vec(K))
            file_tg = File_i(file_sr,IRANK)
            Open (Unit = 10, File=File_tg, status='old', ACTION='read')
            Read(10,*) Seed_vec
            Call Ranset(Seed_vec)
            do NT = 1,LTROT
               do I = 1,Size(Op_V,1)
                  Read(10,*) NSIGMA(I,NT) 
               enddo
            enddo
            close(10)
            Deallocate(Seed_vec)
         else
            If (Irank == 0) then
               Write(6,*) 'No initial configuration'
               OPEN(UNIT=5,FILE='seeds',STATUS='old',ACTION='read',IOSTAT=ierr)
               IF (ierr /= 0) THEN
                  WRITE(*,*) 'unable to open <seeds>',ierr
                  STOP
               END IF
               DO I = Isize-1,1,-1
                  Read (5,*) Seed_in
                  CALL MPI_SEND(Seed_in,1,MPI_INTEGER, I, I+1024,MPI_COMM_WORLD,IERR)
               enddo
               Read(5,*) Seed_in
               CLOSE(5)
            else
               CALL MPI_RECV(Seed_in, 1, MPI_INTEGER,0,  IRANK + 1024,  MPI_COMM_WORLD,STATUS,IERR)
            endif
            Allocate (Seed_vec(1))
            Seed_vec(1) = Seed_in
            Call Ranset(Seed_vec)
            Deallocate(Seed_vec)
            do NT = 1,LTROT
               do I = 1,Size(Op_V,1)
                  X = RANF_WRAP()
                  NSIGMA(I,NT) = 1
                  IF (X.GT.0.5) NSIGMA(I,NT) = -1
               enddo
            enddo
         endif
            
#else   
         INQUIRE (FILE='confin_0', EXIST=lconf)
         If (lconf) Then
            file_tg = "confin_0"
            Call Get_seed_Len(K)
            Allocate(Seed_vec(K))
            Open (Unit = 10, File=File_tg, status='old', ACTION='read')
            Read(10,*) Seed_vec
            Call Ranset(Seed_vec)
            do NT = 1,LTROT
               do I = 1,Size(Op_V,1)
                  Read(10,*) NSIGMA(I,NT) 
               enddo
            enddo
            close(10)
            Deallocate(Seed_vec)
         else
            Write(6,*) 'No initial configuration'
            OPEN(UNIT=5,FILE='seeds',STATUS='old',ACTION='read',IOSTAT=ierr)
            IF (ierr /= 0) THEN
               WRITE(*,*) 'unable to open <seeds>',ierr
               STOP
            END IF
            Read (5,*) Seed_in
            CLOSE(5)
            Allocate(Seed_vec(1))
            Seed_vec(1) = Seed_in
            Call Ranset(Seed_vec)
            Deallocate(Seed_vec)
            do NT = 1,LTROT
               do I = 1,Size(Op_V,1)
                  X = RANF_WRAP()
                  NSIGMA(I,NT) = 1
                  IF (X.GT.0.5) NSIGMA(I,NT) = -1
               enddo
            enddo
         endif
#endif
         
       END SUBROUTINE CONFIN
