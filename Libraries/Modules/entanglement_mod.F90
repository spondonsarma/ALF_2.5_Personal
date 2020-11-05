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

      INTERFACE Calc_Renyi_Ent
        MODULE PROCEDURE Calc_Renyi_Ent_indep, Calc_Renyi_Ent_gen_fl, Calc_Renyi_Ent_gen_all
      END INTERFACE
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
        ! ! (n, n+1), with n=0,2,...
        ! Subroutine Calc_Mutual_Inf(GRC,List_c,Nsites_c,List_f,Nsites_f,Renyi_c,Renyi_f,Renyi_cf)

        !   Implicit none
          
        !   Complex (Kind=8), INTENT(IN)      :: GRC(:,:,:)
        !   Integer, Dimension(:), INTENT(IN) :: List_c, List_f
        !   Integer, INTENT(IN)               :: Nsites_c ,Nsites_f
        !   Complex (Kind=8), INTENT(OUT)   :: Renyi_c, Renyi_f, Renyi_cf

        !   Integer, Dimension(:), Allocatable :: List_cf
        !   Integer          :: I, J, IERR, INFO, Nsites_cf

        !   Nsites_cf=Nsites_c+Nsites_f
          
        !   allocate(List_cf(Nsites_cf))
          
        !   DO I = 1, Nsites_c
        !      List_cf(I) = List_c(I)
        !   END DO
        !   DO I = 1, Nsites_f
        !      List_cf(I+Nsites_c) = List_f(I)
        !   END DO
          
        !   Call Calc_Renyi_Ent(GRC,List_c,Nsites_c,Renyi_c)
        !   Call Calc_Renyi_Ent(GRC,List_f,Nsites_f,Renyi_f)
        !   Call Calc_Renyi_Ent(GRC,List_cf,Nsites_cf,Renyi_cf)
          
        !   deallocate(List_cf)
          
        ! End Subroutine Calc_Mutual_Inf

          Subroutine Calc_Renyi_Ent_indep(GRC,List,Nsites,N_SUN,Renyi)
#ifdef MPI
          Use mpi
#endif
  
          Implicit none
          
          Complex (Kind=8), INTENT(IN)      :: GRC(:,:,:)
          Integer, INTENT(IN)               :: List(:)
          Integer, INTENT(IN)               :: Nsites, N_SUN
          Complex (Kind=8), INTENT(OUT)   :: Renyi

          Complex (Kind=8), Dimension(:,:), Allocatable :: GreenA, GreenA_tmp, IDA
          ! Integer, Dimension(:), Allocatable :: PIVOT
          Complex (Kind=8) :: DET, PRODDET, alpha, beta
          Integer          :: J, IERR, INFO, N_FL, nf, N_FL_half
          Integer         , Dimension(:,:), Allocatable :: List_tmp
          Integer         , Dimension(2)              :: Nsites_tmp,nf_list,N_SUN_tmp
          
          EXTERNAL ZGEMM
          EXTERNAL ZGETRF
          
          N_FL = size(GRC,3)
          N_FL_half = N_FL/2
            
          Renyi=CMPLX(1.d0,0.d0,kind(0.d0))
          alpha=CMPLX(2.d0,0.d0,kind(0.d0))
          beta =CMPLX(1.d0,0.d0,kind(0.d0))

          if (Nsites==0) then
            Renyi=0.d0
            return
          endif
            
#ifdef MPI
          ! Check if entanglement replica group is of size 2 such that the second reny entropy can be calculated
          if(ENT_SIZE==2) then
          
            allocate(List_tmp(NSITES,2))
            Allocate(GreenA(Nsites,2*Nsites),GreenA_tmp(Nsites,2*Nsites),IDA(Nsites,Nsites)) ! new

            DO J = 1, 2
              List_tmp(:,J)=List(:)
              Nsites_tmp(J) = Nsites
              N_SUN_tmp(J) = N_SUN
            enddo

            DO nf=1,N_FL_half
              
              DO J = 1, 2
                nf_list(J) = 2*nf-2+J
              enddo
              call Calc_Renyi_Ent_pair(GRC,List_tmp,Nsites_tmp,nf_list,N_SUN_tmp,PRODDET,GreenA, GreenA_tmp, IDA)
              Renyi = Renyi * PRODDET
              
            Enddo
              
            if (N_FL/=2*N_FL_half) then
            
              call Calc_Renyi_Ent_single(GRC,List_tmp(:,1),Nsites,N_fl,N_SUN,PRODDET,GreenA, GreenA_tmp, IDA)
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
              
          End Subroutine Calc_Renyi_Ent_indep
        
        
        Subroutine Calc_Renyi_Ent_gen_fl(GRC,List,Nsites,N_SUN,Renyi)
#ifdef MPI
          Use mpi
#endif

          Implicit none
          
          Complex (Kind=8), INTENT(IN)      :: GRC(:,:,:)
          Integer, Dimension(:,:), INTENT(IN) :: List ! new
          Integer, INTENT(IN)               :: Nsites(:), N_SUN(:) ! new
          Complex (Kind=8), INTENT(OUT)   :: Renyi

          Complex (Kind=8), Dimension(:,:), Allocatable :: GreenA, GreenA_tmp, IDA
          ! Integer, Dimension(:), Allocatable :: PIVOT
          Complex (Kind=8) :: DET, PRODDET, alpha, beta
          Integer          :: I, J, IERR, INFO, N_FL, nf, N_FL_half, x, dim, dim_eff, nf_eff, start_flav
          Integer         , Dimension(:), Allocatable :: SortedFlavors ! new
          Integer         , Dimension(:,:), Allocatable :: List_tmp
          Integer         , Dimension(2)              :: Nsites_tmp,nf_list,N_SUN_tmp
          
          EXTERNAL ZGEMM
          EXTERNAL ZGETRF
          
          N_FL = size(GRC,3)

          Allocate(SortedFlavors(N_FL)) ! new

          ! insertion sort for small number of elements, SortedFlavors lists flavors in order of size
          start_flav=0
          if (Nsites(1)==0) start_flav = 1
          DO I=2,N_FL
            x = Nsites(I)
            if (Nsites(I)==0) start_flav = start_flav + 1
            J = I-1
            DO while(J >= 1)
              if(Nsites(J) <= x) exit
              SortedFlavors(J+1) = J 
              J = J - 1
            end do
            SortedFlavors(J+1) = I
          END DO
          if(start_flav==N_FL) then
            Renyi=0.0d0
            return
          endif
          N_FL_half = (N_FL-start_flav)/2
          
          Renyi=CMPLX(1.d0,0.d0,kind(0.d0))
          alpha=CMPLX(2.d0,0.d0,kind(0.d0))
          beta =CMPLX(1.d0,0.d0,kind(0.d0))
          
#ifdef MPI
          ! Check if entanglement replica group is of size 2 such that the second reny entropy can be calculated
          if(ENT_SIZE==2) then
          
            !Allocate(GreenA(Nsites,2*Nsites),GreenA_tmp(Nsites,2*Nsites),IDA(Nsites,Nsites))
            dim = Nsites(SortedFlavors(N_FL)) ! new
            allocate(List_tmp(dim,2))
            Allocate(GreenA(dim,2*dim),GreenA_tmp(dim,2*dim),IDA(dim,dim)) ! new

            DO nf=1,N_FL_half

              DO J = 1, 2
                nf_eff = SortedFlavors(start_flav+2*nf-2+J)
                Nsites_tmp(J)=Nsites(nf_eff)
                List_tmp(:,J)=List(:,nf_eff)
                N_sun_tmp(J)=N_SUN(nf_eff)
                nf_list(J)=nf_eff
              enddo
              call Calc_Renyi_Ent_pair(GRC,List_tmp,Nsites_tmp,nf_list,N_SUN_tmp,PRODDET,GreenA, GreenA_tmp, IDA)
              Renyi = Renyi * PRODDET
              
            Enddo
              
            if (N_FL/=2*N_FL_half+start_flav) then

              nf_eff = SortedFlavors(N_fl)
              List_tmp(:,1)=List(:,nf_eff)
            
              call Calc_Renyi_Ent_single(GRC,List_tmp(:,1),Nsites(nf_eff),nf_eff,N_SUN(nf_eff),PRODDET,GreenA, GreenA_tmp, IDA)
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
            
        End Subroutine Calc_Renyi_Ent_gen_fl

        Subroutine Calc_Renyi_Ent_gen_all(GRC,List,Nsites,Renyi)
#ifdef MPI
          Use mpi
#endif

          Implicit none
          
          Complex (Kind=8), INTENT(IN)      :: GRC(:,:,:)
          Integer, Dimension(:,:,:), INTENT(IN) :: List
          Integer, INTENT(IN)               :: Nsites(:,:)
          Complex (Kind=8), INTENT(OUT)   :: Renyi

          Complex (Kind=8), Dimension(:,:), Allocatable :: GreenA, GreenA_tmp, IDA
          ! Integer, Dimension(:), Allocatable :: PIVOT
          Complex (Kind=8) :: DET, PRODDET, alpha, beta
          Integer          :: I, J, IERR, INFO, N_FL, nf, N_FL_half, x, dim, dim_eff, nf_eff, start_flav
          Integer          :: nc, num_nc
          Integer         , Dimension(:), Allocatable :: SortedFlavors,N_SUN_fl,df_list
          Integer         , Dimension(:,:), Allocatable :: List_tmp, eff_ind, eff_ind_inv
          Integer         , Dimension(2)              :: Nsites_tmp,nf_list,N_SUN_tmp
          
          EXTERNAL ZGEMM
          EXTERNAL ZGETRF
          
          N_FL = size(GRC,3)
          num_nc = size(List,3)
          Allocate(SortedFlavors(num_nc*N_FL),N_SUN_fl(num_nc*N_FL),eff_ind(N_FL,num_nc),eff_ind_inv(2,N_FL*num_nc))

          I=0
          do nc=1,num_nc
            do nf=1,N_FL
              I=I+1
              eff_ind(nf,nc)=i
              eff_ind_inv(1,I)=nf
              eff_ind_inv(2,I)=nc
            enddo
          enddo
          N_SUN_fl=1

          ! insertion sort for small number of elements, SortedFlavors lists flavors in order of size
          start_flav=0
          if (Nsites(1,1)==0) start_flav = 1
          ! might have an update in the future to exchange color and flavor loops--optimization
          DO I=2,N_FL*num_nc
            x = Nsites(eff_ind_inv(1,I),eff_ind_inv(2,I))
            if (Nsites(eff_ind_inv(1,I),eff_ind_inv(2,I))==0) start_flav = start_flav + 1 ! there was a bug here
            J = I-1
            DO while(J >= 1)
              if(Nsites(eff_ind_inv(1,J),eff_ind_inv(2,J)) <= x) exit
              SortedFlavors(J+1) = J 
              J = J - 1
            end do
            SortedFlavors(J+1) = I
          END DO
          if(start_flav==N_FL*num_nc) then
            Renyi=0.0d0
            return
          endif

          N_FL_half = (N_FL*num_nc-start_flav)/2
          
          Renyi=CMPLX(1.d0,0.d0,kind(0.d0))
          alpha=CMPLX(2.d0,0.d0,kind(0.d0))
          beta =CMPLX(1.d0,0.d0,kind(0.d0))
            
#ifdef MPI
          ! Check if entanglement replica group is of size 2 such that the second reny entropy can be calculated
          if(ENT_SIZE==2) then
          
            !Allocate(GreenA(Nsites,2*Nsites),GreenA_tmp(Nsites,2*Nsites),IDA(Nsites,Nsites))
            nf=eff_ind_inv(1,SortedFlavors(N_FL*num_nc))
            nc=eff_ind_inv(2,SortedFlavors(N_FL*num_nc))
            dim = Nsites(nf,nc) ! new
            allocate(List_tmp(dim,2))
            Allocate(GreenA(dim,2*dim),GreenA_tmp(dim,2*dim),IDA(dim,dim)) ! new

            DO I=1,N_FL_half

              DO J = 1, 2
                nf=eff_ind_inv(1,SortedFlavors(start_flav+2*I-2+J))
                nc=eff_ind_inv(2,SortedFlavors(start_flav+2*I-2+J))
                Nsites_tmp(J)=Nsites(nf,nc)
                List_tmp(:,J)=List(:,nf,nc)
                N_sun_tmp(J)=1
                nf_list(J)=nf
              enddo
              call Calc_Renyi_Ent_pair(GRC,List_tmp,Nsites_tmp,nf_list,N_SUN_tmp,PRODDET,GreenA, GreenA_tmp, IDA)
              Renyi = Renyi * PRODDET
              
            Enddo
              
            if (N_FL*num_nc/=2*N_FL_half+start_flav) then

              nf=eff_ind_inv(1,SortedFlavors(N_FL*num_nc))
              nc=eff_ind_inv(2,SortedFlavors(N_FL*num_nc))
              List_tmp(:,1)=List(:,nf,nc)
            
              call Calc_Renyi_Ent_single(GRC,List_tmp(:,1),Nsites(nf,nc),nf,1,PRODDET,GreenA, GreenA_tmp, IDA)
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

            
        End Subroutine Calc_Renyi_Ent_gen_all
        
#ifdef MPI
        subroutine Calc_Renyi_Ent_pair(GRC,List,Nsites,nf_list,N_SUN,Renyi,GreenA, GreenA_tmp, IDA)
          Use mpi

          Implicit none
          
          Complex (Kind=8), INTENT(IN)      :: GRC(:,:,:)
          Integer, Dimension(:,:), INTENT(IN) :: List ! new
          Integer, INTENT(IN)               :: Nsites(2), N_SUN(2),nf_list(2) ! new
          Complex (Kind=8), INTENT(OUT), Dimension(:,:), Allocatable :: GreenA, GreenA_tmp, IDA
          Complex (Kind=8), INTENT(OUT)   :: Renyi

          Integer, Dimension(:), Allocatable :: PIVOT
          Complex (Kind=8) :: DET, PRODDET, alpha, beta
          Integer          :: I, J, IERR, INFO, N_FL, nf, N_FL_half, x, dim, dim_eff, nf_eff, start_flav
          Integer         , Dimension(:), Allocatable :: SortedFlavors ! new


          Renyi=CMPLX(1.d0,0.d0,kind(0.d0))
          alpha=CMPLX(2.d0,0.d0,kind(0.d0))
          beta =CMPLX(1.d0,0.d0,kind(0.d0))
          
          dim=size(IDA,1)
          ! We store the reduced Green's function in GreenA_c_tmp (c electrons)
          ! and GreenA_f_tmp (f electrons)
          ! The first Nsites columns contains the spin up sector,
          ! the last Nsites column the spin down sector
          nf_eff = nf_list(1)

          DO J = 1, Nsites(1)
            Do I = 1, Nsites(1)
              GreenA_tmp(I,J) = GRC(List(I,1),List(J,1),nf_eff)
            END DO
          END Do

          nf_eff = nf_list(2)
          DO J = 1, Nsites(2)
            DO I = 1, Nsites(2)
                GreenA_tmp(I,J+dim) = GRC(List(I,2), List(J,2), nf_eff)
            END DO
          END DO
          ! This exchange the last Nsites columns of GreenA_c_tmp between the two replicas
          ! such that GreenA contains the reduced Green's function for two replicas
          ! and a fixed spin sector.
          CALL MPI_ALLTOALL(GreenA_tmp, dim**2, MPI_COMPLEX16, GreenA, dim**2, MPI_COMPLEX16, ENTCOMM, IERR)

          ! Compute Identity - GreenA(replica=1) - GreenA(replica=2) + 2 GreenA(replica=1) * GreenA(replica=2)
          dim_eff = Nsites(1+ENT_RANK)
          IDA(1:dim_eff,1:dim_eff) = - GreenA(1:dim_eff, 1:dim_eff) - GreenA(1:dim_eff, dim+1:dim+dim_eff)
          DO I = 1, dim_eff
              IDA(I,I) = IDA(I,I) + CMPLX(1.d0,0.d0,kind(0.d0))
          END DO
          CALL ZGEMM('n', 'n', dim_eff, dim_eff, dim_eff, alpha, GreenA(1, 1), &
              & dim, GreenA(1, dim+1), dim, beta, IDA, dim)
          ! Compute determinant
          SELECT CASE(dim_eff)
          CASE (1)
            DET = IDA(1,1)
          CASE (2)
            DET = IDA(1,1) * IDA(2,2) - IDA(1,2) * IDA(2,1)
          CASE DEFAULT
            Allocate(PIVOT(dim_eff))
            CALL ZGETRF(dim_eff, dim_eff, IDA, dim, PIVOT, INFO)
            DET = cmplx(1.D0,0.D0,KIND(0.D0))
            DO I = 1, dim_eff
                IF (PIVOT(I).NE.I) THEN
                  DET = -DET * IDA(I,I)
                ELSE
                  DET = DET * IDA(I,I)
                END IF
            ENDDO
            Deallocate(PIVOT)
          END SELECT
          DET=DET**N_SUN(1+ENT_RANK)  !! ATTENTION take care of this!!!
          ! Compute the product of determinants for up and down spin sectors.
          CALL MPI_ALLREDUCE(DET, PRODDET, 1, MPI_COMPLEX16, MPI_PROD, ENTCOMM, IERR)
          ! Now each thread contains in PRODDET the full determinant, as obtained by
          ! a pair of replicas.
          Renyi = Renyi * PRODDET
        end subroutine Calc_Renyi_Ent_pair

        subroutine Calc_Renyi_Ent_single(GRC,List,Nsites,nf_eff,N_SUN,Renyi,GreenA, GreenA_tmp, IDA)
          Use mpi
          
          Implicit none
          
          Complex (Kind=8), INTENT(IN)      :: GRC(:,:,:)
          Integer, Dimension(:), INTENT(IN) :: List ! new
          Integer, INTENT(IN)               :: Nsites, N_SUN,nf_eff ! new
          Complex (Kind=8), INTENT(OUT), Dimension(:,:), Allocatable :: GreenA, GreenA_tmp, IDA
          Complex (Kind=8), INTENT(OUT)   :: Renyi

          Integer, Dimension(:), Allocatable :: PIVOT
          Complex (Kind=8) :: DET, PRODDET, alpha, beta
          Integer          :: I, J, IERR, INFO, N_FL, nf, N_FL_half, x, dim, dim_eff, start_flav
          Integer         , Dimension(:), Allocatable :: SortedFlavors ! new


          Renyi=CMPLX(1.d0,0.d0,kind(0.d0))
          alpha=CMPLX(2.d0,0.d0,kind(0.d0))
          beta =CMPLX(1.d0,0.d0,kind(0.d0))
          
          dim=size(IDA,1)
          ! We store the reduced Green's function in GreenA_c_tmp (c electrons)
          ! and GreenA_f_tmp (f electrons)
          ! The first Nsites columns contains the last flavour sector,
          ! the last Nsites column are old (Nf>2) or random, but they won't be used

          DO J = 1, Nsites
            DO I = 1, Nsites
                GreenA_tmp(I,J) = GRC(List(I), List(J), nf_eff)
            END DO
          END DO
          
          ! This exchange the last Nsites columns of GreenA_c_tmp between the two replicas
          ! such that GreenA contains the reduced Green's function for two replicas
          ! and a fixed spin sector.
          CALL MPI_ALLTOALL(GreenA_tmp, dim**2, MPI_COMPLEX16, GreenA, dim**2, MPI_COMPLEX16, ENTCOMM, IERR)
          
          DET = cmplx(1.D0,0.D0,KIND(0.D0))
          if(ENT_RANK==0) then
            dim_eff = NSites
            ! Compute Identity - GreenA(replica=1) - GreenA(replica=2) + 2 GreenA(replica=1) * GreenA(replica=2)
            IDA = - GreenA(1:dim_eff, 1:dim_eff) - GreenA(1:dim_eff, dim+1:dim+dim_eff)
            DO I = 1, dim_eff
                IDA(I,I) = IDA(I,I) + CMPLX(1.d0,0.d0,kind(0.d0))
            END DO
            CALL ZGEMM('n', 'n', dim_eff, dim_eff, dim_eff, alpha, GreenA(1, 1), &
                & dim, GreenA(1, dim+1), dim, beta, IDA, dim)
            ! Compute determinant
            SELECT CASE(dim_eff)
            CASE (1)
              DET = IDA(1,1)
            CASE (2)
              DET = IDA(1,1) * IDA(2,2) - IDA(1,2) * IDA(2,1)
            CASE DEFAULT
              Allocate(PIVOT(dim_eff))
              CALL ZGETRF(dim_eff, dim_eff, IDA, dim, PIVOT, INFO)
              DO I = 1, dim_eff
                  IF (PIVOT(I).NE.I) THEN
                    DET = -DET * IDA(I,I)
                  ELSE
                    DET = DET * IDA(I,I)
                  END IF
              ENDDO
              Deallocate(PIVOT)
            END SELECT
            Det=Det**N_SUN
          endif
          
          ! Compute the product of determinants for up and down spin sectors.
          CALL MPI_ALLREDUCE(DET, PRODDET, 1, MPI_COMPLEX16, MPI_PROD, ENTCOMM, IERR)
          ! Now each thread contains in PRODDET the full determinant, as obtained by
          ! a pair of replicas.
          Renyi = Renyi * PRODDET
        end subroutine Calc_Renyi_Ent_single
#endif
        
      end Module entanglement
