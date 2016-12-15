Module Random_Wrap

   contains

     Subroutine Get_seed_Len(K)
       Implicit none
       Integer :: K
       CALL RANDOM_SEED (SIZE=K) 
     end Subroutine Get_seed_Len

     
     Subroutine Ranset(Iseed_vec)
       Implicit none
       Integer, Dimension(:) :: Iseed_vec

       Integer :: K, N, i, Iseed
       Integer, allocatable :: Seed_start(:)
       Real (Kind=8) :: X

       N = size(Iseed_vec) 
       CALL RANDOM_SEED (SIZE=K)           
       Allocate         (SEED_start(K) ) 
       ! Setup SEED_start
       Iseed = Iseed_vec(1)
       If (Iseed == 0 )  then 
          Iseed = 8752143
          N = 0
       endif
       do i = 1,K
          if (i <= N) then 
             SEED_Start(i) = Iseed_vec(i)
          else
             X = lcg(Iseed)
             SEED_Start(i) = Iseed
          endif
       enddo
       CALL RANDOM_SEED (PUT = SEED_start(1:K)) 
       !Write(6,*) 'Starting seeds ', SEED_Start

     end Subroutine Ranset
       
     Subroutine Ranget(Iseed_vec)
       Implicit none
       Integer, Dimension(:) :: Iseed_vec

       Integer :: K, N, i, Iseed
       Integer, allocatable :: Seed_end(:)
       Real (Kind=8) :: X

       N = size(Iseed_vec) 
       CALL RANDOM_SEED (SIZE=K)           
       Allocate         (SEED_end(K) ) 
       CALL RANDOM_SEED (GET = SEED_end(1:K)) 
       ! Setup SEED_start
       Iseed = Iseed_vec(1)
       do i = 1,N
          if (i <= K) then 
             Iseed_vec(i)  =  SEED_end(i)
          else
             X = lcg(Iseed)
             Iseed_vec(i) = Iseed
          endif
       enddo
       !Write(6,*) 'End seeds ', SEED_end

     end Subroutine Ranget

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function returns a real Pseudo-Random Number using a Linear congrential
!> Random number generator. The range of the returned values is [0.0, 1.0).
!> Note that the value 1.0 will not be returned.
!
!> @param[in] seed An integer to seed the LCG.
!-------------------------------------------------------------------- 
     real (Kind=Kind(0.D0)) function lcg(seed)
       implicit none
       integer :: seed
       integer(8) :: res, norm
       
       res = seed ! convert type
       res = 62089911*res + 4349
       norm = 2147483648_8 !specify 8 byte integer
       lcg = DBLE(MODULO(res, 2147483647))/DBLE(norm)
     end function lcg

     real (Kind=8)  function  ranf(iq)
       implicit none
       integer, optional ::  iq
       Real (Kind=8) :: X
       Call Random_Number(X)
       ranf = X
     end function ranf


     real (Kind=8)  function  ranf_wrap(iq)
       implicit none
       integer, optional ::  iq
       Real (Kind=8) :: X
       Call Random_Number(X)
       ranf_wrap = X
     end function ranf_wrap
     

      real (kind=8)  function  rang(iq)

        ! Random variable according to the distribution:  exp(-x**2/2)/(sqrt(2*3.1415927))
      
        integer iq
        real (Kind=8) ::  pi, ranmod, theta
      
        PI = 3.1415926536D0
        RANMOD = SQRT(-2.D0 * LOG(RANF(iq)))
        THETA  = 2.D0 * PI * RANF(iq)
        rang = RANMOD * COS(THETA)
        
      end function rang
   
   end Module Random_Wrap
