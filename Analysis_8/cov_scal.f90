      Program Cov_vec

        
        Use ERRORS
        Implicit none

        REAL    (KIND=8), DIMENSION(:,:), ALLOCATABLE :: OBS
        REAL    (KIND=8), DIMENSION(:),   ALLOCATABLE :: EN, SIGN
        REAL    (KIND=8) :: XM, XERR

        Complex (Kind=8) Z1,Z2,Z3,Z4,Z5
        Complex (Kind=8), Allocatable  :: Tmp(:)
        Integer :: Nobs 
        Integer :: NST, NS, NS1, NS2, NSTEP, NC, NP,  Nbins, NP_EFF, ISEED, I, IOBS
        Integer :: N,N1, NBIN

        ! Count the number of bins

        Open (Unit=10, File="Var_scal", status="unknown")
        Read(10,*) NOBS
        allocate (Tmp(NOBS) )
        rewind(10)
        nbins = 0
        do
            read(10,*,End=10) N,  (Tmp(I), I=1,size(Tmp,1))
            nbins = nbins + 1
         enddo
10       continue
         Write(6,*) "# of bins: ", Nbins
         Close(10) 
         
         NP = NBINS
         
         ALLOCATE(OBS(NP,NOBS))

         !Open (Unit=25, File="statdat1", status="unknown") 
         !read(25,*) NST, NS1, NS2, NSTEP
         !Close(25)
         NST = 2; NS1 = 1; NS2 = 1; NSTEP = 1
         OPEN (UNIT=20, FILE='Var_scal', STATUS='old')
         NC = 0
         DO N = 1,NP
            IF (N.GE.NST) THEN
               NC = NC + 1
               READ(20,*) N1, (Tmp(I), I=1,size(Tmp,1)) 
               OBS(NC,:) = dble(Tmp(:)) 
            ELSE
               READ(20,*) N1, (Tmp(I), I=1,size(Tmp,1)) 
            ENDIF
         ENDDO
         CLOSE(20)
2100     FORMAT(I6,2X,F16.8)
         
         OPEN (UNIT=21, FILE='Var_scalJ', STATUS='unknown')
         WRITE(21,*) 'Effective number of bins, and bins: ', NC, NP
         NP_EFF = NC
         ALLOCATE (EN(NP_EFF), SIGN(NP_EFF))
         DO IOBS = 1,NOBS
            WRITE(21,*)
            DO I = 1,NP_EFF
               EN  (I) = OBS(I,IOBS)
               SIGN(I) = OBS(I,NOBS)
            ENDDO
            DO NBIN = NS1, NS2, NSTEP
               if (NBIN.gt.0) then
                  IF (IOBS.EQ.NOBS  ) then 
                     CALL ERRCALCJ(EN,XM,XERR,NBIN)
                  else
                     CALL ERRCALCJ(EN,SIGN,XM,XERR,NBIN)
                  endif
                  WRITE(21,2001) IOBS, XM,  XERR
                  ! Test
                  ! NBOOT = 40
                  ! CALL BOOTSTRAP( EN,XM_BS,XERR_BS,NBOOT,ISEED)
                  ! WRITE(21,2001) IOBS, XM_BS,  XERR_BS
                  ! IF (IOBS == 4) Write(22,"(F14.7,2x,F14.7)")  XM/dble(L*L), XERR/dble(L*L)
               endif
            ENDDO
         ENDDO
         CLOSE(21)
2001     FORMAT('OBS : ', I4,4x,F12.6,2X, F12.6)
         
         DEALLOCATE (EN,SIGN,OBS)
         
       END Program Cov_vec
       
