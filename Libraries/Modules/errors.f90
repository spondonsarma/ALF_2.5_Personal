!  Copyright (C) 2016 The ALF project
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

     MODULE ERRORS
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Collection of routines for error analysis 
!
!--------------------------------------------------------------------

       Use MyMats
       Use Random_Wrap

       INTERFACE ERRCALC
          MODULE PROCEDURE ERRCALC, ERRCALC_C
       END INTERFACE
       INTERFACE ERRCALCJ
          MODULE PROCEDURE ERRCALC_J,    ERRCALC_J_REBIN,    ERRCALC_JS, ERRCALC_JS_REBIN, &
               &           ERRCALC_J_C,  ERRCALC_J_C_REBIN,  ERRCALC_JS_C, ERRCALC_JS_C_REBIN
       END INTERFACE
       INTERFACE COV
          MODULE PROCEDURE COVJ, COVJS, COVJS_C, COVJS_C_REBIN
       END INTERFACE
       INTERFACE COV_ERR
          MODULE PROCEDURE COV_ERR
       END INTERFACE
       INTERFACE INTERGRATE_F
          MODULE PROCEDURE INTER_F
       END INTERFACE
       INTERFACE INTERGRATE
          MODULE PROCEDURE INTER_QMC
       END INTERFACE
       INTERFACE FIT
          MODULE PROCEDURE FIT
       END INTERFACE
       INTERFACE AUTO_COR
          MODULE PROCEDURE  AUTO_COR
       END INTERFACE
       INTERFACE Bootstrap
          MODULE PROCEDURE  Bootstrap
       END INTERFACE
       INTERFACE Bootstrap_fluc
          MODULE PROCEDURE  BootstrapC_fluc
       END INTERFACE

       

       CONTAINS 
!***********
         SUBROUTINE ERRCALC(EN,XM,XERR)
!	  Calculates error on the input vector EN.  Just the standard deviation.

           IMPLICIT NONE
           REAL (Kind=Kind(0.d0)), DIMENSION(:) ::  EN
           REAL (Kind=Kind(0.d0))               ::  XM, XERR, XSQ
           INTEGER                     ::  NP, NT
 
           NP = SIZE(EN)
           
           XM  = 0.D0
           DO NT = 1,NP
              XM  = XM  + EN(NT)
           ENDDO
           XM    = XM /DBLE(NP)
           XSQ = 0.D0
           DO NT = 1,NP
              XSQ = XSQ + (EN(NT)-XM)**2
           ENDDO
           XSQ   = XSQ/DBLE(NP)
           XERR  = XSQ/DBLE(NP)
           IF (XERR.GT.0.D0) THEN
              XERR = SQRT(XERR)
           ELSE
              XERR = 0.D0 
           ENDIF
           
           RETURN
         END SUBROUTINE ERRCALC


         SUBROUTINE ERRCALC_C(EN,ZM,ZERR)
!	  Calculates error on the input vector EN.  Just the standard deviation.

           IMPLICIT NONE
           Complex (Kind=Kind(0.d0)), DIMENSION(:) ::  EN
           Complex (Kind=Kind(0.d0))               ::  ZM, ZERR
           INTEGER                     ::  NP, NT
 
           ! Local 
           Real (Kind=Kind(0.d0)), dimension(:), allocatable :: Rhelp
           real (Kind=Kind(0.d0)) :: XM, XERR

           NP = SIZE(EN)
           Allocate (Rhelp(NP))
           
           do nt = 1,np
              Rhelp(nt) = dble(en(nt))
           enddo
           call errcalc(Rhelp, xm, xerr)
           zm   =  cmplx(xm  , 0.d0, kind(0.D0))
           Zerr =  cmplx(xerr, 0.d0, kind(0.D0))

           do nt = 1,np
              Rhelp(nt) = aimag(en(nt))
           enddo
           call errcalc(Rhelp, xm, xerr)
           zm   =  zm   + cmplx( 0.d0, xm, kind(0.D0)   )
           Zerr =  Zerr + cmplx( 0.d0, xerr, kind(0.D0) )

           RETURN
         END SUBROUTINE ERRCALC_C

         SUBROUTINE ERRCALC_J(EN,XM,XERR)
!	   Calculates jacknife error on the input vector EN.  Mean and  variance.
!          The input are the bins.
           
           IMPLICIT NONE

           REAL (Kind=Kind(0.d0)), DIMENSION(:) ::  EN
           REAL (Kind=Kind(0.d0))               ::  XM, XERR, X, Xhelp
           REAL (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE ::  EN1
           INTEGER     :: NP, N, N1

           NP = SIZE(EN)
           ALLOCATE (EN1(NP))
           
           ! Build the jackknife averages and send to errcalc.

           Xhelp = 0.D0
           DO N1 = 1,NP
              Xhelp = Xhelp + EN(N1)
           ENDDO
              
           DO N = 1,NP
              X = Xhelp - EN(N)
              EN1(N) = X / DBLE(NP -1)
           ENDDO
           CALL ERRCALC(EN1,XM,XERR)
           XERR = XERR*DBLE(NP)
           DEALLOCATE  ( EN1 )

           RETURN
         END SUBROUTINE ERRCALC_J


         SUBROUTINE ERRCALC_J_C(EN,ZM,ZERR)
!	   Calculates jacknife error on the input vector EN.  Mean and  variance.
!          The input are the bins.
           
           IMPLICIT NONE

           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:) ::  EN
           COMPLEX (Kind=Kind(0.d0))               ::  ZM, ZERR, Z, Zhelp
           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE ::  EN1
           INTEGER     :: NP, N, N1

           NP = SIZE(EN)
           ALLOCATE (EN1(NP))
           
           ! Build the jackknife averages and send to errcalc.

           Zhelp = CMPLX(0.D0, 0.D0, kind(0.D0))
           DO N1 = 1,NP
              Zhelp = Zhelp + EN(N1)
           ENDDO
              
           DO N = 1,NP
              Z =  Zhelp - EN(N)
              EN1(N) = Z / DBLE(NP -1)
           ENDDO
           CALL ERRCALC(EN1,ZM,ZERR)
           ZERR = ZERR*DBLE(NP)
           DEALLOCATE  ( EN1 )

           RETURN
         END SUBROUTINE ERRCALC_J_C

!************
         SUBROUTINE ERRCALC_J_C_REBIN(EN,ZM,ZERR,NREBIN)
!	   Calculates jacknife error on the input vector EN.  Mean and  variance.
!          The input are the bins.
           
           IMPLICIT NONE

           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:) ::  EN
           COMPLEX (Kind=Kind(0.d0))               ::  ZM, ZERR, Z
           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE ::  EN1
           INTEGER     :: NP, N, NP1, NREBIN, NC, NB

           NP = SIZE(EN)
           NP1 = NP/NREBIN
           ALLOCATE (EN1(NP1))
           ! Rebin
           NC = 0
           DO N = 1,NP1
              Z = CMPLX(0.D0, 0.D0, kind(0.D0))
              DO NB = 1,NREBIN
                 NC = NC + 1
                 Z = Z + EN(NC)
              ENDDO
              Z = Z/DBLE(NREBIN)
              EN1(N) = Z
           ENDDO
           CALL ERRCALC_J_C(EN1,ZM,ZERR)
           
           DEALLOCATE(EN1)
           
         END SUBROUTINE ERRCALC_J_C_REBIN

!******************
         SUBROUTINE ERRCALC_J_REBIN(EN,XM,XERR,NREBIN)
!	   Calculates jacknife error on the input vector EN with rebinning.  Mean and  variance.
!          The input are the bins.

           IMPLICIT NONE

           REAL (Kind=Kind(0.d0)), DIMENSION(:) ::  EN
           REAL (Kind=Kind(0.d0))               ::  XM, XERR, X
           REAL (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE ::  EN1
           INTEGER :: NREBIN, NC, N, NB, NP1, NP
 
           NP = SIZE(EN)
           NP1 = NP/NREBIN
           ALLOCATE (EN1(NP1))
           
           ! Rebin
           NC = 0
           DO N = 1,NP1
              X = 0.D0
              DO NB = 1,NREBIN
                 NC = NC + 1
                 X = X + EN(NC)
              ENDDO
              X = X/DBLE(NREBIN)
              EN1(N) = X
           ENDDO
           CALL ERRCALC_J(EN1,XM,XERR)

           DEALLOCATE(EN1)
           RETURN
         END SUBROUTINE ERRCALC_J_REBIN

!**********
         SUBROUTINE ERRCALC_JS(EN,SI,XM,XERR)
!	   Calculates error on the input vector EN.  Just the variance.
!          The input are the bins

           IMPLICIT NONE

           REAL (Kind=Kind(0.d0)), DIMENSION(:) ::  EN, SI
           REAL (Kind=Kind(0.d0))               ::  XM, XERR, X,XS
           REAL (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE :: EN1
           INTEGER                     ::  N, N1, NP, NP1

           NP = SIZE(EN)
           NP1= SIZE(SI)
           IF (NP1.NE.NP) THEN
              WRITE(6,*) 'Error in Errcalc_JS'
              STOP
           ENDIF
           ALLOCATE (EN1(NP))
           
           ! Build the jackknife averages and send to errcalc

           DO N = 1,NP
              X  = 0.D0
              XS = 0.D0
              DO N1 = 1,NP
                 IF (N1.NE.N)  X  = X  + EN(N1)
                 IF (N1.NE.N)  XS = XS + SI(N1)
              ENDDO
              EN1(N) = X / XS
           ENDDO
           CALL ERRCALC(EN1,XM,XERR)
           XERR = XERR*DBLE(NP)
           DEALLOCATE  ( EN1 )

           RETURN
         END SUBROUTINE ERRCALC_JS

!**********
         SUBROUTINE ERRCALC_JS_C(EN,SI,XM,XERR)
!	   Calculates error on the input vector EN.  Just the variance.
!          The input are the bins

           IMPLICIT NONE

           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:) ::  EN, SI
           COMPLEX (Kind=Kind(0.d0))               ::  XM, XERR, X,XS
           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE :: EN1
           INTEGER                     ::  N, N1, NP, NP1

           NP = SIZE(EN)
           NP1= SIZE(SI)
           IF (NP1.NE.NP) THEN
              WRITE(6,*) 'Error in Errcalc_JS'
              STOP
           ENDIF
           ALLOCATE (EN1(NP))
           
           ! Build the jackknife averages and send to errcalc

           DO N = 1,NP
              X  = CMPLX(0.D0, 0.D0, kind(0.D0))
              XS = CMPLX(0.D0, 0.D0, kind(0.D0))
              DO N1 = 1,NP
                 IF (N1.NE.N)  X  = X  + EN(N1)
                 IF (N1.NE.N)  XS = XS + SI(N1)
              ENDDO
              EN1(N) = X / XS
           ENDDO
           CALL ERRCALC(EN1,XM,XERR)
           XERR = XERR*DBLE(NP)
           DEALLOCATE  ( EN1 )

           RETURN
         END SUBROUTINE ERRCALC_JS_C



!********
         SUBROUTINE ERRCALC_JS_REBIN(EN,SI,XM,XERR,NREBIN)
!	   Calculates jacknife error on the input vector EN with rebinning.  Mean and  variance.
!          The input are the bins.
           
           IMPLICIT NONE

           REAL (Kind=Kind(0.d0)), DIMENSION(:) ::  EN, SI
           REAL (Kind=Kind(0.d0))               ::  XM, XERR, X, Y
           REAL (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE ::  EN1, SI1
           INTEGER :: NREBIN, NC, N, NB, NP, NP1
 
           NP = SIZE(EN)
           NP1 = NP/NREBIN
           ALLOCATE (EN1(NP1))
           ALLOCATE (SI1(NP1))
           
           ! Rebin
           NC = 0
           DO N = 1,NP1
              X = 0.D0; Y = 0.D0
              DO NB = 1,NREBIN
                 NC = NC + 1
                 X = X + EN(NC)
                 Y = Y + SI(NC)
              ENDDO
              X = X/DBLE(NREBIN)
              Y = Y/DBLE(NREBIN)
              EN1(N) = X
              SI1(N) = Y
           ENDDO
           CALL ERRCALC_JS(EN1,SI1,XM,XERR)

           DEALLOCATE (EN1,SI1)

           RETURN
         END SUBROUTINE ERRCALC_JS_REBIN

!******************
         SUBROUTINE ERRCALC_JS_C_REBIN(EN,SI,XM,XERR,NREBIN)
!	   Calculates jacknife error on the input vector EN with rebinning.  Mean and  variance.
!          The input are the bins.
           
           IMPLICIT NONE

           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:) ::  EN, SI
           COMPLEX (Kind=Kind(0.d0))               ::  XM, XERR, X, Y
           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE ::  EN1, SI1
           INTEGER :: NREBIN, NC, N, NB, NP, NP1
 
           NP = SIZE(EN)
           NP1 = NP/NREBIN
           ALLOCATE (EN1(NP1))
           ALLOCATE (SI1(NP1))
           
           ! Rebin
           NC = 0
           DO N = 1,NP1
              X = cmplx(0.D0,0.d0,kind(0.d0)); Y = cmplx(0.D0,0.D0,kind(0.d0))
              DO NB = 1,NREBIN
                 NC = NC + 1
                 X = X + EN(NC)
                 Y = Y + SI(NC)
              ENDDO
              X = X/DBLE(NREBIN)
              Y = Y/DBLE(NREBIN)
              EN1(N) = X
              SI1(N) = Y
           ENDDO
           CALL ERRCALC_JS_C(EN1,SI1,XM,XERR)

           DEALLOCATE (EN1,SI1)

           RETURN
         END SUBROUTINE ERRCALC_JS_C_REBIN

!******************

         SUBROUTINE INTER_QMC(GR, SIGN1, DTAU, RES, ERR) 
           
           IMPLICIT NONE
           ! Given GR(Times, Bins)  and Sign1(Bins) calculates the integral and error
           ! The sign is the same for all Times.
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:) ::  GR
           REAL (Kind=Kind(0.d0)), DIMENSION(:)   ::   SIGN1

           !Local
           REAL (Kind=Kind(0.d0)), DIMENSION(:  ), ALLOCATABLE  ::  HLP
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:), ALLOCATABLE  ::  HLP1
           REAL (Kind=Kind(0.d0))                 ::  X, Y, Err, Res, DTAU
           INTEGER :: NT, NB, NB1, NTDM, NDATA
           
           NTDM  = SIZE(GR,1)
           NDATA = SIZE(GR,2)


           ALLOCATE( HLP(NDATA), HLP1(NTDM,NDATA) )
           DO NT = 1,NTDM
              DO NB= 1, NDATA
                 X = 0.D0
                 Y = 0.D0
                 DO NB1 = 1,NDATA
                    IF (NB1.NE.NB) THEN
                       X = X + GR(NT,NB1)
                       Y = Y + SIGN1(NB1)
                    ENDIF
                 ENDDO
                 HLP1(NT,NB) = X/Y
              ENDDO
           ENDDO
           
           DO NB = 1,NDATA
              X = 0.D0
              DO NT = 1,NTDM-1
                 X = X + (HLP1(NT,NB) + HLP1(NT+1,NB))*0.5D0
              ENDDO
              HLP (NB   )  = X * DTAU
           ENDDO
           
           CALL ERRCALC(HLP, RES, ERR) 
           ERR = ERR*DBLE(NDATA)

           DEALLOCATE( HLP, HLP1 )

           RETURN
         END SUBROUTINE INTER_QMC
         
!******************
         REAL (Kind=Kind(0.d0)) FUNCTION INTER_F(A,B,N,F)
           ! integrates the function F from A to B  using N points.
           
           IMPLICIT NONE

           INTEGER  :: N, I
           REAL (Kind=Kind(0.d0)) ::  A, B, X, X1
           REAL (Kind=Kind(0.d0)), EXTERNAL :: F

           REAL (Kind=Kind(0.d0)) ::  DEL

           DEL = (B-A)/DBLE(N)
           INTER_F = 0.D0
           DO I = 0, N-1
              X  = A + DBLE(I  )*DEL
              X1 = A + DBLE(I+1)*DEL 
              INTER_F = INTER_F + ( F(X) + F(X1) )*0.5D0  
           ENDDO
           INTER_F = INTER_F*DEL
         END FUNCTION INTER_F

!****************** Least square fits:
         SUBROUTINE FIT(XDATA,FDATA,ERROR,ARES,CHSQ,F)

           IMPLICIT NONE

           REAL (Kind=Kind(0.d0)), DIMENSION(:) ::  XDATA, FDATA, ERROR, ARES
           REAL (Kind=Kind(0.d0))               ::  CHSQ,  X
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:),  ALLOCATABLE :: A, U,V,VINV,V1
           REAL (Kind=Kind(0.d0)), DIMENSION(:  ),  ALLOCATABLE :: B,D
           REAL (Kind=Kind(0.d0)), EXTERNAL  :: F
           INTEGER                     ::  NDATA, NBASIS, I, M, M1, NCON, N

           NDATA = SIZE(XDATA)
           NBASIS= SIZE(ARES)

           !WRITE(6,*) 'NDATA, NBASIS: ',NDATA, NBASIS 
           ALLOCATE (A(NDATA,NBASIS))
           ALLOCATE (U(NDATA,NBASIS))
           ALLOCATE (D(NBASIS))
           ALLOCATE (V   (NBASIS,NBASIS))
           ALLOCATE (V1  (NBASIS,NBASIS))
           ALLOCATE (VINV(NBASIS,NBASIS))
           ALLOCATE (B(NDATA))
           
           A = 0.D0
           U = 0.D0
           D = 0.D0
           V = 0.D0
           VINV = 0.D0
           V1 = 0.D0
           B = 0.D0
           NCON = 1
           DO M = 1,NBASIS
              DO I = 1,NDATA
                 A(I,M) = F(M,XDATA(I))/ERROR(I)
              ENDDO
           ENDDO
           DO I = 1,NDATA
              B(I) = FDATA(I)/ERROR(I)
           ENDDO
           !write(6,*) A
           CALL UDV(A,U,D,V,NCON)
           DO M = 1,NBASIS
              DO I = 1,NBASIS
                 V1(I,M) = V(M,I)
              ENDDO
           ENDDO
           X = 0.D0
           CALL INV(V1,VINV,X)

           DO M1 = 1,NBASIS
              X = 0.D0
              DO M = 1,NBASIS
                 DO I = 1,NDATA
                    X = X + B(I)*U(I,M)*VINV(M,M1)/D(M)
                 ENDDO
              ENDDO
              ARES(M1) = X
           ENDDO

           CHSQ = 0.D0
           DO N = 1,NDATA
              X = 0.D0
              DO M = 1,NBASIS
                 X = X + ARES(M)*F(M,XDATA(N))
              ENDDO
              CHSQ = CHSQ + (FDATA(N) - X)**2/ERROR(N)**2
           ENDDO
           CHSQ = CHSQ/DBLE(NDATA)
           
           DEALLOCATE (A)
           DEALLOCATE (U)
           DEALLOCATE (D)
           DEALLOCATE (V)
           DEALLOCATE (V1)
           DEALLOCATE (VINV)
           DEALLOCATE (B)
           
         END SUBROUTINE FIT
	
         SUBROUTINE COVJ(GR, XCOV, XMEAN) 
           
           IMPLICIT NONE
           !Given GR(Times, Bins)  calculates the mean and the covariance.
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:) ::  GR, XCOV
           REAL (Kind=Kind(0.d0)), DIMENSION(:)   ::  XMEAN

           !Local
           REAL (Kind=Kind(0.d0)), DIMENSION(:  ), ALLOCATABLE  ::  HLP
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:), ALLOCATABLE  ::  HLP1
           REAL (Kind=Kind(0.d0))                 ::  X, XM, XERR
           INTEGER :: NT, NT1, NB, NB1, NTDM, NDATA
           
           NTDM  = SIZE(GR,1)
           NDATA = SIZE(GR,2)

           IF ( (SIZE(XCOV,1).NE.SIZE(XCOV,2) ) .OR. (SIZE(XCOV,1).NE.NTDM) ) THEN
              WRITE(6,*) 'Error in COV'
              STOP
           ENDIF

           ALLOCATE( HLP(NDATA), HLP1(NTDM,NDATA) )
           DO NT = 1,NTDM
              DO NB= 1, NDATA
                 X = 0.0
                 DO NB1 = 1,NDATA
                    IF (NB1.NE.NB) THEN
                       X = X + GR(NT,NB1)
                    ENDIF
                 ENDDO
                 HLP1(NT,NB) = X/DBLE(NDATA-1)
                 HLP (NB   ) = X/DBLE(NDATA-1)
              ENDDO
              CALL ERRCALC(HLP,XM ,XERR)
              XMEAN(NT) = XM
           ENDDO
	
           
           DO NT = 1,NTDM
              DO NT1= 1,NTDM
                 X = 0.0
                 DO NB = 1,NDATA
                    X = X +  HLP1(NT,NB)*HLP1(NT1,NB)
                 ENDDO
                 X = X/DBLE(NDATA)
                 XCOV(NT,NT1)  = ( X - XMEAN(NT)*XMEAN(NT1) )*DBLE(NDATA)
              ENDDO
           ENDDO
           

           DEALLOCATE( HLP, HLP1 )

           RETURN
         END SUBROUTINE COVJ


         SUBROUTINE COVJS(GR, SIGN1, XCOV, XMEAN) 
           
           IMPLICIT NONE
           ! Given GR(Times, Bins)  and Sign1(Bins) calculates the mean and the covariance.
           ! The sign is the same for all Times.
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:) ::  GR, XCOV
           REAL (Kind=Kind(0.d0)), DIMENSION(:)   ::  XMEAN, SIGN1

           !Local
           REAL (Kind=Kind(0.d0)), DIMENSION(:  ), ALLOCATABLE  ::  HLP
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:), ALLOCATABLE  ::  HLP1
           REAL (Kind=Kind(0.d0))                 ::  X, XM, XERR, Y
           INTEGER :: NT, NT1, NB, NB1, NTDM, NDATA
           
           NTDM  = SIZE(GR,1)
           NDATA = SIZE(GR,2)

           IF ( (SIZE(XCOV,1).NE.SIZE(XCOV,2) ) .OR. (SIZE(XCOV,1).NE.NTDM) ) THEN
              WRITE(6,*) 'Error in COV'
              STOP
           ENDIF

           ALLOCATE( HLP(NDATA), HLP1(NTDM,NDATA) )
           DO NT = 1,NTDM
              DO NB= 1, NDATA
                 X = 0.D0
                 Y = 0.D0
                 DO NB1 = 1,NDATA
                    IF (NB1.NE.NB) THEN
                       X = X + GR(NT,NB1)
                       Y = Y + SIGN1(NB1)
                    ENDIF
                 ENDDO
                 HLP1(NT,NB) = X/Y
                 HLP (NB   ) = X/Y
              ENDDO
              CALL ERRCALC(HLP,XM ,XERR)
              XMEAN(NT) = XM
           ENDDO
	
           
           DO NT = 1,NTDM
              DO NT1= 1,NTDM
                 X = 0.0
                 DO NB = 1,NDATA
                    X = X +  HLP1(NT,NB)*HLP1(NT1,NB)
                 ENDDO
                 X = X/DBLE(NDATA)
                 XCOV(NT,NT1)  = ( X - XMEAN(NT)*XMEAN(NT1) )*DBLE(NDATA)
              ENDDO
           ENDDO
           

           DEALLOCATE( HLP, HLP1 )

           RETURN
         END SUBROUTINE COVJS




         SUBROUTINE COVJS_C(GR, SIGN1, XCOV, XMEAN) 
           
           IMPLICIT NONE
           ! Given GR(Times, Bins)  and Sign1(Bins) calculates the mean and the covariance.
           ! The sign is the same for all Times. 
           Complex (Kind=Kind(0.d0)), DIMENSION(:,:) ::  GR, XCOV
           Complex (Kind=Kind(0.d0)), DIMENSION(:)   ::  XMEAN
           Real    (Kind=Kind(0.d0)), DIMENSION(:)   ::  SIGN1


           !Local
           REAL (Kind=Kind(0.d0)), DIMENSION(:  ), ALLOCATABLE  ::  HLP, XMEAN_R
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:), ALLOCATABLE  ::  HLP1
           REAL (Kind=Kind(0.d0))                 ::  X, XM, XERR, Y
           INTEGER :: NT, NT1, NB, NB1, NTDM, NDATA, Nth
           COMPLEX (Kind=Kind(0.d0)) :: Z

           NTDM  = SIZE(GR,1)
           NDATA = SIZE(GR,2)


           !Write(6,*) 'Errors.f90 ', NTDM, NDATA
           IF ( (SIZE(XCOV,1).NE.SIZE(XCOV,2) ) .OR. (SIZE(XCOV,1).NE.NTDM) ) THEN
              WRITE(6,*) 'Error in COV'
              STOP
           ENDIF


           
           ALLOCATE( HLP(NDATA), HLP1(NTDM,NDATA), XMEAN_R(NTDM) )
           XMEAN = CMPLX(0.d0, 0.d0, kind(0.D0))
           XCOV  = CMPLX(0.d0, 0.d0, kind(0.D0))
           
           DO NTH = 1,2
              Z = CMPLX(1.D0, 0.D0, kind(0.D0))
              IF (NTH .EQ. 2 ) Z = CMPLX( 0.D0, -1.D0, kind(0.D0))
              DO NT = 1,NTDM
                 DO NB= 1, NDATA
                    X = 0.D0
                    Y = 0.D0
                    DO NB1 = 1,NDATA
                       IF (NB1.NE.NB) THEN
                          X = X + DBLE ( Z*GR(NT,NB1) ) 
                          Y = Y + SIGN1(NB1)
                       ENDIF
                    ENDDO
                    HLP1(NT,NB) = X/Y
                    HLP (NB   ) = X/Y
                 ENDDO
                 CALL ERRCALC(HLP,XM ,XERR)
                 XMEAN(NT) = XMEAN(NT) + CONJG(Z)*XM
                 XMEAN_R(NT) = XM
                 !if (Nth.eq.2) write(6,*) XM
              ENDDO
              
              
              DO NT = 1,NTDM
                 DO NT1= 1,NTDM
                    X = 0.0
                    DO NB = 1,NDATA
                       X = X +  HLP1(NT,NB)*HLP1(NT1,NB)
                    ENDDO
                    X = X/DBLE(NDATA)
                    XCOV(NT,NT1)  = XCOV(NT,NT1) + CONJG(Z)* ( X - XMEAN_R(NT)*XMEAN_R(NT1) )*DBLE(NDATA)
                 ENDDO
              ENDDO
           ENDDO

           DEALLOCATE( HLP, HLP1, XMEAN_R )
           
           RETURN
         END SUBROUTINE COVJS_C
         

!========================

         SUBROUTINE COVJS_C_REBIN(GR, SIGN2, XCOV, XMEAN,NREBIN) 
           
           IMPLICIT NONE
           ! Given GR(Times, Bins)  and Sign1(Bins) calculates the mean and the covariance.
           ! The sign is the same for all Times. 
           Complex (Kind=Kind(0.d0)), DIMENSION(:,:) ::  GR, XCOV
           Complex (Kind=Kind(0.d0)), DIMENSION(:)   ::  XMEAN
           Real    (Kind=Kind(0.d0)), DIMENSION(:)   ::  SIGN2
           INTEGER :: NREBIN


           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:,:), ALLOCATABLE ::  GR1
           REAL    (Kind=Kind(0.d0)), DIMENSION(:), ALLOCATABLE ::  SIGN1
           
           INTEGER :: NTDM, NDATA, NDATA1, N, NB, NC, NT
           REAL    (Kind=Kind(0.d0)) :: X
           COMPLEX (Kind=Kind(0.d0)) :: Z
           
           NTDM  = SIZE(GR,1)
           NDATA = SIZE(GR,2)

           NDATA1 = NDATA/NREBIN
           ALLOCATE ( GR1(NTDM,NDATA1), SIGN1(NDATA1) )
           
           SIGN1 = 0.d0
           GR1   = CMPLX(0.d0,0.d0,kind(0.d0))
           ! Rebin
           NC = 0
           DO N = 1,NDATA1
              DO NB = 1,NREBIN
                 NC = NC + 1
                 SIGN1(N) = SIGN1(N) + SIGN2(NC)
                 DO NT = 1,NTDM
                    GR1(NT,N)  = GR1(NT,N) + GR(NT,NC)
                 ENDDO
              ENDDO
           ENDDO
           SIGN1 = SIGN1/DBLE (NREBIN)
           GR1   = GR1  /CMPLX(DBLE(NREBIN),0.d0,KIND(0.d0))

           CALL COVJS_C(GR1, SIGN1, XCOV, XMEAN) 
           
           DEALLOCATE ( GR1, SIGN1 )

         END SUBROUTINE COVJS_C_REBIN


         Subroutine COV_ERR(XMEAN, XCOV, ISEED) 
           !  Given Mean and Cov, diagonalizes the COV and produces a new data set within 
           !  the errorbars
           
           Implicit None
           ! Parameters
           REAL (Kind=Kind(0.d0)), DIMENSION(:,:) :: XCOV
           REAL (Kind=Kind(0.d0)), DIMENSION(:)   :: XMEAN 
           
           Integer :: ntau, I, M, ISeed, ISEED_VEC(1)
           Real (Kind=Kind(0.d0)) :: X

           Real (Kind=Kind(0.d0)), Dimension(:,:),  allocatable ::  UC
           Real (Kind=Kind(0.d0)), Dimension(:),    allocatable ::  XMEAN_1, SIG_1

           ntau = size(Xmean,1) 
           Allocate (UC(ntau,ntau), XMEAN_1(ntau), SIG_1(ntau) ) 

           ISEED_VEC(1) = ISEED
           CALL RANSET(ISEED_VEC)

           CALL DIAG(XCOV,UC,SIG_1)
           
           DO I = 1,NTAU
              X = 0.D0
              DO M = 1,NTAU
                 X  = X + UC(M,I)* XMEAN(M)
              ENDDO
              XMEAN_1(I) = X
           ENDDO
           DO I = 1,NTAU
              IF (SIG_1(I).LT.0.d0) Then 
                  write(6,*) 'Error in Cov_err', SIG_1(I)
              Endif
              XMEAN_1(I) = XMEAN_1(I) + SQRT(ABS(SIG_1(I)))*RANG_WRAP()
           ENDDO
           DO I = 1,NTAU
              X = 0.D0
              DO M = 1,NTAU
                 X  = X + UC(I,M)*XMEAN_1(M)
              ENDDO
              XMEAN(I) = X
           ENDDO

           CALL RANGET(ISEED_VEC)
           ISEED  = ISEED_VEC(1)

           Deallocate (UC, XMEAN_1, SIG_1) 


         END Subroutine COV_ERR
         
         SUBROUTINE  AUTO_COR(DATA,RES)

           Implicit none
           
           REAL (Kind=Kind(0.d0)),  DIMENSION(:)  :: DATA,RES

           !Local 
           Integer  :: nb, nt, ntau, nt1
           Real (Kind=Kind(0.d0)) :: X1, X2, X3
          
           nb = SIZE(DATA)
           nt = SIZE(RES) 
           if (nb.lt.nt) then
              write(6,*) 'Error in autocor'
              stop
           end if
           
           DO ntau = 1,  nt
              X1 = 0.0
              X2 = 0.0
              X3 = 0.0
              DO nt1 = 1, nb - ntau
                 X3 = X3 + DATA(nt1)
              ENDDO
              X3 = X3 / dble(nb - ntau)

              DO nt1 = 1, nb - ntau
                 X1 = X1 + (DATA(nt1)-x3)*(DATA(nt1 + ntau)-x3) 
                 X2 = X2 + (DATA(nt1)-x3)*(DATA(nt1)-x3)
              ENDDO
              X1 = X1 / dble(nb - ntau)
              X2 = X2 / dble(nb - ntau)
              
              Res(ntau)  = X1/X2
              
           ENDDO
           
         END SUBROUTINE AUTO_COR

         SUBROUTINE BOOTSTRAPC_FLUC(A,B,AB,NBOOT,ISEED,ZM,ZERR)
           !!!  COMPUTES <AB> - <A><B>
           IMPLICIT NONE
           COMPLEX (Kind=Kind(0.d0)), DIMENSION(:), INTENT(IN) :: A,B,AB
           INTEGER, INTENT(IN)    :: NBOOT
           INTEGER, INTENT(INOUT) :: ISEED
           COMPLEX (Kind=Kind(0.d0)), INTENT(OUT) :: ZM,ZERR
           

           !Local
           INTEGER :: NP, NB, I, J, ISEED_VEC(1)
           COMPLEX (Kind=Kind(0.d0)) :: Z,  Z1,Z2,Z12
           

           ISEED_VEC(1) = ISEED
           CALL RANSET(Iseed_vec)
           NP = SIZE(A,1) 
           ZM   = CMPLX(0.d0,0.d0,Kind=Kind(0.d0))
           ZERR = CMPLX(0.d0,0.d0,Kind=Kind(0.d0))
           DO NB = 1, NBOOT
              Z1  = cmplx(0.d0,0.d0,Kind=Kind(0.d0))
              Z2  = cmplx(0.d0,0.d0,Kind=Kind(0.d0))
              Z12 = cmplx(0.d0,0.d0,Kind=Kind(0.d0))
              DO I = 1,NP
                 J = NINT( DBLE(NP)* RANF_WRAP() + 0.5 )
                 IF (J == 0) J = 1
                 IF (J > NP) J = NP
                 Z1 = Z1  + A(J)
                 Z2 = Z2  + B(J)
                 Z12 =Z12 + AB(J)
              ENDDO
              Z1 = Z1 /CMPLX(DBLE(NP),0.d0,Kind=Kind(0.d0))
              Z2 = Z2 /CMPLX(DBLE(NP),0.d0,Kind=Kind(0.d0)) 
              Z12 =Z12/CMPLX(DBLE(NP),0.d0,Kind=Kind(0.d0))

              Z    = Z12 - Z1*Z2
              ZM   = ZM   + Z
              ZERR = ZERR + Z*Z
           ENDDO
           ZM   = ZM  /CMPLX(DBLE(NBOOT),0.d0,Kind=Kind(0.d0))
           ZERR = ZERR/CMPLX(DBLE(NBOOT),0.d0,Kind=Kind(0.d0))
           
           Z = ZERR -  ZM*ZM
           ZERR = SQRT(Z)
           
           CALL RANGET(Iseed_vec)
           ISEED = ISEED_VEC(1)


         END SUBROUTINE BOOTSTRAPC_FLUC

         SUBROUTINE BOOTSTRAP(EN,XM,XERR,NBOOT,ISEED)

           IMPLICIT NONE
           REAL (Kind=Kind(0.d0)), DIMENSION(:) ::  EN
           REAL (Kind=Kind(0.d0))               ::  XM, XERR,  X
           INTEGER                     ::  NP, NT, NBOOT, NB, I, ISEED
           
           ! Local 
           INTEGER :: ISEED_VEC(1)
 
           NP = SIZE(EN)
           ISEED_VEC(1) = ISEED
           CALL RANSET(Iseed_vec)
           
           
           ! Build the Bootstrap samples
           
           XM   = 0.D0
           XERR = 0.D0
           DO NB = 1,NBOOT
              X = 0.D0
              DO NT = 1, NP
                 I = NINT( DBLE(NP)* RANF_WRAP() + 0.5 )
                 IF (I.EQ.0 .OR. I.GT.NP ) THEN
                    WRITE(6,*) 'ERROR IN BOOTSTRAP'
                    STOP
                 ENDIF
                 X = X + EN(I)
              ENDDO
              X = X/DBLE(NP)
              XM   = XM + X
              XERR = XERR + X*X
           ENDDO

           XM   = XM  /DBLE(NBOOT)
           XERR = XERR/DBLE(NBOOT)

           X = XERR - XM*XM
           XERR = 0.d0
           IF (X.GT.0.d0) XERR = SQRT(X)

           CALL RANGET(Iseed_vec)
           ISEED = ISEED_VEC(1) 

         END SUBROUTINE BOOTSTRAP

       END MODULE ERRORS


