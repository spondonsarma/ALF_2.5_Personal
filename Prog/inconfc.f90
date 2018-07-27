!  Copyright (C) 2016 - 2018 The ALF project
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

        SUBROUTINE CONFIN

!--------------------------------------------------------------------
!
!> @brief 
!> Reads in the field configuration and seeds if present so as to 
!> pursue a run. If  the configuration is not present  the 
!> routine will generate one based on the seeds read in from the file
!> seeds
!
!--------------------------------------------------------------------
         USE HAMILTONIAN
#ifdef MPI
         USE mpi
#endif
         IMPLICIT NONE

         ! LOCAL
         INTEGER        :: I, IERR, SEED_IN, K, NT
         INTEGER, DIMENSION(:), ALLOCATABLE :: SEED_VEC
         LOGICAL ::   LCONF 
         CHARACTER (LEN=64) :: FILE_SR, FILE_TG, FILE_seeds, FILE_info, File1

#ifdef MPI
         INTEGER        :: STATUS(MPI_STATUS_SIZE), irank_g, isize_g, igroup, ISIZE, IRANK
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
         call MPI_Comm_rank(Group_Comm, irank_g, ierr)
         call MPI_Comm_size(Group_Comm, isize_g, ierr)
         igroup           = irank/isize_g
#endif

         ALLOCATE (NSIGMA(SIZE(OP_V,1),LTROT))
         
#if defined(MPI)

#if defined(TEMPERING) 
            write(FILE1,'(A,I0,A)') "Temp_",igroup,"/confin_0"
#else
            File1 = "confin_0"
#endif
            INQUIRE (FILE=File1, EXIST=LCONF)
            IF (LCONF) THEN
#if defined(TEMPERING) 
               write(FILE_TG,'(A,I0,A,I0)') "Temp_",igroup,"/confin_",irank_g
#else
               write(FILE_TG,'(A,I0)') "confin_",irank_g
#endif
               CALL GET_SEED_LEN(K)
               ALLOCATE(SEED_VEC(K))
               OPEN (UNIT = 10, FILE=FILE_TG, STATUS='OLD', ACTION='READ')
               READ(10,*) SEED_VEC
               CALL RANSET(SEED_VEC)
               DO NT = 1,LTROT
                  DO I = 1,SIZE(OP_V,1)
                     READ(10,*) NSIGMA(I,NT) 
                  ENDDO
               ENDDO
               CLOSE(10)
               DEALLOCATE(SEED_VEC)
            ELSE
               IF (IRANK == 0) THEN
                  WRITE(6,*) 'No initial configuration'
                  OPEN(UNIT=5,FILE='seeds',STATUS='OLD',ACTION='READ',IOSTAT=IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'UNABLE TO OPEN <seeds>',IERR
                     STOP
                  END IF
                  DO I = ISIZE-1,1,-1
                     READ (5,*) SEED_IN
                     CALL MPI_SEND(SEED_IN,1,MPI_INTEGER, I, I+1024, MPI_COMM_WORLD,IERR)
                  ENDDO
                  READ(5,*) SEED_IN
                  CLOSE(5)
               ELSE
                  CALL MPI_RECV(SEED_IN, 1, MPI_INTEGER,0,  IRANK + 1024,  MPI_COMM_WORLD,STATUS,IERR)
               ENDIF
               ALLOCATE (SEED_VEC(1))
               SEED_VEC(1) = SEED_IN
               CALL RANSET(SEED_VEC)
               DEALLOCATE (SEED_VEC)
               Call Hamiltonian_set_random_nsigma
               if (irank_g == 0) then
#if defined(TEMPERING) 
                  write(FILE_info,'(A,I0,A)') "Temp_",igroup,"/info"
#else
                  FILE_info="info"
#endif
                  Open (Unit = 50,file=FILE_info,status="unknown",position="append")
                  WRITE(50,*) 'No initial configuration, Seed_in', SEED_IN
                  Close(50)
               endif
         ENDIF
            
#else   
         FILE_TG = "confin_0"
         INQUIRE (FILE=FILE_TG, EXIST=LCONF)
         IF (LCONF) THEN
            CALL GET_SEED_LEN(K)
            ALLOCATE(SEED_VEC(K))
            OPEN (UNIT = 10, FILE=FILE_TG, STATUS='OLD', ACTION='READ')
            READ(10,*) SEED_VEC
            CALL RANSET(SEED_VEC)
            DO NT = 1,LTROT
               DO I = 1,SIZE(OP_V,1)
                  READ(10,*) NSIGMA(I,NT) 
               ENDDO
            ENDDO
            CLOSE(10)
            DEALLOCATE(SEED_VEC)
         ELSE
            FILE_seeds="seeds"
            OPEN(UNIT=5,FILE=FILE_seeds,STATUS='OLD',ACTION='READ',IOSTAT=IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'UNABLE TO OPEN <seeds>',IERR
               STOP
            END IF
            READ (5,*) SEED_IN
            CLOSE(5)
            FILE_info="info"
            Open (Unit = 50,file=FILE_info,status="unknown",position="append")
            WRITE(50,*) 'No initial configuration, Seed_in', SEED_IN
            Close(50)

            ALLOCATE(SEED_VEC(1))
            SEED_VEC(1) = SEED_IN
            CALL RANSET (SEED_VEC)
            DEALLOCATE  (SEED_VEC)
            Call Hamiltonian_set_random_nsigma
            
         ENDIF
#endif
         
       END SUBROUTINE CONFIN
