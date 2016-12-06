      Program Cov_vec

        
        Use ERRORS
        Implicit none

        REAL    (KIND=8), DIMENSION(:,:), ALLOCATABLE :: OBS
        REAL    (KIND=8), DIMENSION(:),   ALLOCATABLE :: EN, SIGN
        REAL    (KIND=8) :: XM, XERR, X

        Complex (Kind=8) Z1,Z2,Z3,Z4,Z5
        Complex (Kind=8), Allocatable  :: Tmp(:)
        Integer :: Nobs 
        Integer :: NST, NS, NS1, NS2, NSTEP, NC, NP,  Nbins, NP_EFF, ISEED, I, IOBS
        Integer :: N,N1, NBIN

        Integer :: n_skip, N_rebin, N_Cov, ierr
        NAMELIST /VAR_errors/   n_skip, N_rebin, N_Cov


         OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'unable to open <parameters>',ierr
            STOP
         END IF
         READ(5,NML=VAR_errors)
         CLOSE(5)


        ! Count the number of bins
        Open (Unit=10, File="Var_scal", status="unknown")
        Read(10,*) NOBS
        allocate (Tmp(NOBS) )
        rewind(10)
        nbins = 0
        do
            read(10,*,End=10) N,  (Tmp(I), I=1,size(Tmp,1)-1), X
            Tmp(NOBS) = cmplx(X,0.d0,kind(0.d0))
            nbins = nbins + 1
         enddo
10       continue
         Write(6,*) "# of bins: ", Nbins
         Close(10) 
         
         NP = NBINS
         
         ALLOCATE(OBS(NP,NOBS))

         NST = N_skip; NS1 = N_Rebin; NS2 = N_Rebin; NSTEP = 1
         OPEN (UNIT=20, FILE='Var_scal', STATUS='old')
         NC = 0
         DO N = 1,NP
            IF (N.GE.NST) THEN
               NC = NC + 1
               READ(20,*) N1, (Tmp(I), I=1,size(Tmp,1)-1),X 
               Tmp(NOBS) = cmplx(X,0.d0,kind(0.d0))
               OBS(NC,:) = dble(Tmp(:)) 
            ELSE
               READ(20,*) N1, (Tmp(I), I=1,size(Tmp,1)-1),X 
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
               EN  (I) = Real(OBS(I,IOBS), kind(0.d0))
               SIGN(I) = Real(OBS(I,NOBS), kind(0.d0))
            ENDDO
            DO NBIN = NS1, NS2, NSTEP
               if (NBIN.gt.0) then
                  IF (IOBS.EQ.NOBS  ) then 
                     CALL ERRCALCJ(EN,XM,XERR,NBIN)
                  else
                     CALL ERRCALCJ(EN,SIGN,XM,XERR,NBIN)
                  endif
                  WRITE(21,2001) IOBS, XM,  XERR
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
       
