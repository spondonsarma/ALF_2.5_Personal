!  Copyright (C) 2018 The ALF project
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
 
     Module entanglement

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This module generates one and two dimensional Bravais lattices.
!
!--------------------------------------------------------------------

      ! Used for MPI
      INTEGER :: ENTCOMM, ENT_RANK, ENT_SIZE=0, Norm, group

      Contains
!========================================================================

        Subroutine Init_Entanglement_replicas(Group_Comm)
#ifdef MPI
          Use mpi
#endif
          Implicit none
          Integer, INTENT(IN)               :: Group_Comm
          
          Integer                           :: ISIZE, IRANK, IERR
#ifdef MPI
          ! Create subgroups of two replicas each
          CALL MPI_COMM_SIZE(Group_Comm,ISIZE,IERR)
          CALL MPI_COMM_RANK(Group_Comm,IRANK,IERR)
          CALL MPI_COMM_SPLIT(Group_Comm, IRANK / 2, IRANK, ENTCOMM, IERR)
          CALL MPI_COMM_RANK(ENTCOMM,ENT_RANK,IERR)
          CALL MPI_COMM_SIZE(ENTCOMM,ENT_SIZE,IERR)
          Norm=ISIZE/2 ! number of pairs
          norm=2*norm ! but each task of pair contributes
          group=Group_Comm
#endif
          
        end Subroutine Init_Entanglement_replicas
          
!========================================================================
        ! Calculation of the Renyi entanglement entropy
        ! The algorithm works only for an MPI program
        ! We partition the nodes into groups of 2 replicas:
        ! (n, n+1), with n=0,2,...
        Subroutine Calc_Mutual_Inf(GRC,Phase,Ntau,List_c,Nsites_c,List_f,Nsites_f,Renyi_c,Renyi_f,Renyi_cf)

          Implicit none
          
          Complex (Kind=8), INTENT(IN)      :: GRC(:,:,:)
          Complex (Kind=8), Intent(IN)      :: PHASE
          Integer, INTENT(IN)               :: Ntau
          Integer, Dimension(:), INTENT(IN) :: List_c, List_f
          Integer, INTENT(IN)               :: Nsites_c ,Nsites_f
          Complex (Kind=8), INTENT(OUT)   :: Renyi_c, Renyi_f, Renyi_cf

          Integer, Dimension(:), Allocatable :: List_cf
          Integer          :: I, J, IERR, INFO, Nsites_cf

          Nsites_cf=Nsites_c+Nsites_f
          
          allocate(List_cf(Nsites_cf))
          
          DO I = 1, Nsites_c
             List_cf(I) = List_c(I)
          END DO
          DO I = 1, Nsites_f
             List_cf(I+Nsites_c) = List_f(I)
          END DO
          
          Call Calc_Renyi_Ent(GRC,Phase,Ntau,List_c,Nsites_c,Renyi_c)
          Call Calc_Renyi_Ent(GRC,Phase,Ntau,List_f,Nsites_f,Renyi_f)
          Call Calc_Renyi_Ent(GRC,Phase,Ntau,List_cf,Nsites_cf,Renyi_cf)
          
          deallocate(List_cf)
          
        End Subroutine Calc_Mutual_Inf
        
        
        Subroutine Calc_Renyi_Ent(GRC,Phase,Ntau,List,Nsites,Renyi)
#ifdef MPI
          Use mpi
#endif

          Implicit none
          
          Complex (Kind=8), INTENT(IN)      :: GRC(:,:,:)
          Complex (Kind=8), Intent(IN)      :: PHASE
          Integer, INTENT(IN)               :: Ntau
          Integer, Dimension(:), INTENT(IN) :: List
          Integer, INTENT(IN)               :: Nsites
          Complex (Kind=8), INTENT(OUT)   :: Renyi

          Complex (Kind=8), Dimension(:,:), Allocatable :: GreenA, GreenA_tmp, IDA
          Integer, Dimension(:), Allocatable :: PIVOT
          Complex (Kind=8) :: DET, PRODDET
          Integer          :: I, J, IERR, INFO, N_FL, nf, N_FL_half
          
          EXTERNAL ZGEMM
          EXTERNAL ZGETRF
          
          N_FL = size(GRC,3)
          N_FL_half = N_FL/2
          
          Renyi=CMPLX(1.d0,0.d0,kind(0.d0))
          
#ifdef MPI
          ! Check if entanglement replica group is of size 2 such that the second reny entropy can be calculated
          if(ENT_SIZE==2) then
          
            Allocate(GreenA(Nsites,2*Nsites),GreenA_tmp(Nsites,2*Nsites),IDA(Nsites,Nsites))

            DO nf=1,N_FL_half
              ! We store the reduced Green's function in GreenA_c_tmp (c electrons)
              ! and GreenA_f_tmp (f electrons)
              ! The first Nsites columns contains the spin up sector,
              ! the last Nsites column the spin down sector
              DO J = 1, Nsites
                DO I = 1, Nsites
                    GreenA_tmp(I,J) = GRC(List(I), List(J), 2*nf-1)
                END DO
              END DO
              DO J = 1, Nsites
                DO I = 1, Nsites
                    GreenA_tmp(I,J+Nsites) = GRC(List(I), List(J), 2*nf)
                END DO
              END DO
              ! This exchange the last Nsites columns of GreenA_c_tmp between the two replicas
              ! such that GreenA contains the reduced Green's function for two replicas
              ! and a fixed spin sector.
              CALL MPI_ALLTOALL(GreenA_tmp, Nsites**2, MPI_COMPLEX16, GreenA, Nsites**2, MPI_COMPLEX16, ENTCOMM, IERR)

              ! Compute Identity - GreenA(replica=1) - GreenA(replica=2) + 2 GreenA(replica=1) * GreenA(replica=2)
              IDA = - GreenA(1:Nsites, 1:Nsites) - GreenA(1:Nsites, Nsites+1:2*Nsites)
              DO I = 1, Nsites
                  IDA(I,I) = IDA(I,I) + CMPLX(1.d0,0.d0,kind(0.d0))
              END DO
              CALL ZGEMM('n', 'n', Nsites, Nsites, Nsites, CMPLX(2.D0,0.D0,KIND(0.D0)), GreenA(1:Nsites, 1:Nsites), &
                  & Nsites, GreenA(1:Nsites, Nsites+1:2*Nsites), Nsites, CMPLX(1.D0,0.D0,KIND(0.D0)), IDA, Nsites)
              ! Compute determinant
              SELECT CASE(Nsites)
              CASE (1)
                DET = IDA(1,1)
              CASE (2)
                DET = IDA(1,1) * IDA(2,2) - IDA(1,2) * IDA(2,1)
              CASE DEFAULT
                Allocate(PIVOT(Nsites))
                CALL ZGETRF(Nsites, Nsites, IDA, Nsites, PIVOT, INFO)
                DET = cmplx(1.D0,0.D0,KIND(0.D0))
                DO I = 1, Nsites
                    IF (PIVOT(I).NE.I) THEN
                      DET = -DET * IDA(I,I)
                    ELSE
                      DET = DET * IDA(I,I)
                    END IF
                ENDDO
                Deallocate(PIVOT)
              END SELECT
              ! Compute the product of determinants for up and down spin sectors.
              CALL MPI_ALLREDUCE(DET, PRODDET, 1, MPI_COMPLEX16, MPI_PROD, ENTCOMM, IERR)
              ! Now each thread contains in PRODDET the full determinant, as obtained by
              ! a pair of replicas.
              Renyi = Renyi * PRODDET
              
            Enddo
              
            if (N_FL/=2*N_FL_half) then
            
              ! We store the reduced Green's function in GreenA_c_tmp (c electrons)
              ! and GreenA_f_tmp (f electrons)
              ! The first Nsites columns contains the last flavour sector,
              ! the last Nsites column are old (Nf>2) or random, but they won't be used
              DO J = 1, Nsites
                DO I = 1, Nsites
                    GreenA_tmp(I,J) = GRC(List(I), List(J), N_FL)
                END DO
              END DO
              
              ! This exchange the last Nsites columns of GreenA_c_tmp between the two replicas
              ! such that GreenA contains the reduced Green's function for two replicas
              ! and a fixed spin sector.
              CALL MPI_ALLTOALL(GreenA_tmp, Nsites**2, MPI_COMPLEX16, GreenA, Nsites**2, MPI_COMPLEX16, ENTCOMM, IERR)
              
              DET = cmplx(1.D0,0.D0,KIND(0.D0))
              if(ENT_RANK==0) then

                ! Compute Identity - GreenA(replica=1) - GreenA(replica=2) + 2 GreenA(replica=1) * GreenA(replica=2)
                IDA = - GreenA(1:Nsites, 1:Nsites) - GreenA(1:Nsites, Nsites+1:2*Nsites)
                DO I = 1, Nsites
                    IDA(I,I) = IDA(I,I) + CMPLX(1.d0,0.d0,kind(0.d0))
                END DO
                CALL ZGEMM('n', 'n', Nsites, Nsites, Nsites, CMPLX(2.D0,0.D0,KIND(0.D0)), GreenA(1:Nsites, 1:Nsites), &
                    & Nsites, GreenA(1:Nsites, Nsites+1:2*Nsites), Nsites, CMPLX(1.D0,0.D0,KIND(0.D0)), IDA, Nsites)
                ! Compute determinant
                SELECT CASE(Nsites)
                CASE (1)
                  DET = IDA(1,1)
                CASE (2)
                  DET = IDA(1,1) * IDA(2,2) - IDA(1,2) * IDA(2,1)
                CASE DEFAULT
                  Allocate(PIVOT(Nsites))
                  CALL ZGETRF(Nsites, Nsites, IDA, Nsites, PIVOT, INFO)
                  DO I = 1, Nsites
                      IF (PIVOT(I).NE.I) THEN
                        DET = -DET * IDA(I,I)
                      ELSE
                        DET = DET * IDA(I,I)
                      END IF
                  ENDDO
                  Deallocate(PIVOT)
                END SELECT
              endif
                
              ! Compute the product of determinants for up and down spin sectors.
              CALL MPI_ALLREDUCE(DET, PRODDET, 1, MPI_COMPLEX16, MPI_PROD, ENTCOMM, IERR)
              ! Now each thread contains in PRODDET the full determinant, as obtained by
              ! a pair of replicas.
              
              Renyi = Renyi * PRODDET
            
            endif
            
            Deallocate(GreenA,GreenA_tmp,IDA)
            
          else
            ! if there had been an odd number of task in tempering group / world, set renyi to 0
            Renyi=CMPLX(0.d0,0.d0,kind(0.d0))
          endif
          
          ! average over all pairs of replicas, the single task contributes nothing even so it takes part in the call
          CALL MPI_ALLREDUCE(Renyi, PRODDET, 1, MPI_COMPLEX16, MPI_SUM, group, IERR)
          Renyi=PRODDET/dble(norm)
          ! At this point, each task of the temepering group / world returns the same averaged value of the pairs, including the possible "free"/ unpaired one.
          ! This mechanisms leads to some syncronization, but I (Johannes) am lacking a better way to treat odd number of tasks.
#endif
            
        End Subroutine Calc_Renyi_Ent
        
      end Module entanglement
