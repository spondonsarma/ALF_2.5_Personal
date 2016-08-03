
    MODULE MyMats

       INTERFACE MMULT
          !C = A*B MMULT(C, A, B)
          MODULE PROCEDURE MMULT_R, MMULT_C
       END INTERFACE
       INTERFACE INITD
          MODULE PROCEDURE INITD_R, INITD_C
       END INTERFACE
       INTERFACE COMPARE
          MODULE PROCEDURE COMPARE_R, COMPARE_C
       END INTERFACE
       INTERFACE DET
          MODULE PROCEDURE DET_C
       END INTERFACE DET
       INTERFACE INV
          MODULE PROCEDURE INV_R0, INV_R_Variable, INV_R_VARIABLE_1, INV_R1, INV_R2, INV_C, INV_C1, &
               &        INV_C_Variable  
       END INTERFACE
       INTERFACE UDV
          MODULE PROCEDURE UDV1_R, UDV_C
       END INTERFACE
       INTERFACE QR
          MODULE PROCEDURE QR_C
       END INTERFACE QR
       INTERFACE SVD
          MODULE PROCEDURE SVD_C
       END INTERFACE SVD
       INTERFACE DIAG
          MODULE PROCEDURE DIAG_R, DIAG_I
       END INTERFACE
       INTERFACE DIAG_GEN
          MODULE PROCEDURE DIAG_GEN
       END INTERFACE DIAG_GEN
       INTERFACE SECONDS
          MODULE PROCEDURE SECONDS
       END INTERFACE
     CONTAINS

!*************
       SUBROUTINE DIAG_GEN(Z_MAT,U,W,LR,ICON)
         IMPLICIT NONE
         COMPLEX   (KIND=8), INTENT(IN), DIMENSION(:,:) :: Z_MAT
         CHARACTER (LEN=1),  INTENT(IN)  :: LR
         COMPLEX   (KIND=8), INTENT(INOUT), DIMENSION(:,:) :: U
         COMPLEX   (KIND=8), INTENT(INOUT), DIMENSION(:) :: W
         INTEGER :: ICON

         !!!! Uses Lapack !!!
         ! LR = L  then         U*A   = W*U   Left  eigenvectors
         ! LR = R  then         A*U   = W*U   Right eigenvectors
         

         !  Local space
         INTEGER :: N, LDA, LDVL, LDVR, INFO, LWORK, I, J, M
         CHARACTER (LEN=1) ::   JOBVL, JOBVR
         COMPLEX (KIND=8), ALLOCATABLE, DIMENSION(:,:) :: A, VL, VR
         REAL (KIND=8)   , ALLOCATABLE, DIMENSION(:) :: RWORK
         COMPLEX (KIND=8), ALLOCATABLE, DIMENSION(:) :: WORK
         
         REAL    (KIND=8) :: XMAX, X
         COMPLEX (KIND=8) :: Z
         
         N = SIZE(Z_MAT,1)
         ALLOCATE(A(N,N))
         A = Z_MAT
         LDA = N

         JOBVR  = "N"
         JOBVL  = "N"
         LDVL = 1
         LDVR = 1
         IF (LR =="L") THEN 
            JOBVL ="V"
            LDVL  = N
         ELSEIF (LR =="R") THEN
            JOBVR ="V"
            LDVR = N
         ELSE
            WRITE(6,*) 'Error in DIAG_GEN' 
            STOP
         ENDIF
         ALLOCATE(VL(LDVL,N),  VR(LDVR,N) )
         LWORK = 2*N
         ALLOCATE (WORK(LWORK), RWORK(LWORK) )

         CALL ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, &
              &      WORK, LWORK, RWORK, INFO )
         
         IF (LR=="R")  THEN
            DO I = 1,N
               DO J = 1,N
                  U(I,J) = VR(I,J)
               ENDDO
            ENDDO
         ELSE
            DO I = 1,N
               DO J = 1,N 
                  U(I,J) = CONJG(VL(J,I))
               ENDDO
            ENDDO
         ENDIF
         
         IF (ICON == 1 ) THEN 
            !Test 
            XMAX = 0.d0
            DO I = 1,N
               DO J = 1,N
                  IF (LR=="R")  THEN
                     Z = cmplx(0.d0,0.d0,kind=8)
                     DO M = 1,N
                        Z = Z + Z_MAT(I,M)*U(M,J)
                     ENDDO
                     Z = Z - W(I)*U(I,J)
                     X = SQRT( DBLE(  Z*CONJG(Z) ) )
                  ENDIF
                  IF (LR=="L")  THEN
                     Z = cmplx(0.d0,0.d0,kind=8)
                     DO M = 1,N
                        Z = Z + U(I,M)*Z_MAT(M,J) 
                     ENDDO
                     Z = Z - W(I)*U(I,J)
                     X = SQRT( DBLE(  Z*CONJG(Z) ) )
                  ENDIF
                  IF ( X > XMAX ) XMAX = X
               ENDDO
            ENDDO
            WRITE(6,*) 'Testing Diag_GEN :', XMAX
            !End Test
         ENDIF

         DEALLOCATE(VL, VR)
         DEALLOCATE(WORK, RWORK)
         DEALLOCATE(A)
         

       END SUBROUTINE DIAG_GEN
!*************
       SUBROUTINE MMULT_R(C, A, B)
         IMPLICIT NONE
         REAL (KIND=8), DIMENSION(:,:) :: A,B,C
         REAL (KIND=8) :: X, ALP, BET
         INTEGER I,J, K, N, M, P, LDA, LDB, LDC
         N = SIZE(A,1) ! Rows in A
         M = SIZE(A,2) ! Columns in A
         P = SIZE(B,2) ! Columns in B
         LDA = N; LDB = SIZE(B,1); LDC = SIZE(C,1)

         ALP = 1.D0
         BET = 0.D0


         CALL DGEMM('n','n',N,P,M,ALP,A,LDA,B,LDB,BET,C,LDC)






! WRITE(6,*) 'In real', N,M,P
! DO I = 1,N
! DO J = 1,P
! X = 0.D0
! DO K = 1,M
! X = X + A(I,K)*B(K,J)
! ENDDO
! C(I,J) = X
! ENDDO
! ENDDO
       END SUBROUTINE MMULT_R

       SUBROUTINE MMULT_C(C, A, B)
         IMPLICIT NONE
         COMPLEX (KIND=8), DIMENSION(:,:) :: A,B,C
         COMPLEX (KIND=8) :: ALP, BET
         INTEGER I,J, K, N, M, P, LDA, LDB, LDC

         N = SIZE(A,1)
         M = SIZE(A,2)
         P = SIZE(B,2)
         LDA = N; LDB = SIZE(B,1); LDC = SIZE(C,1)

         ALP = DCMPLX(1.D0,0.D0)
         BET = DCMPLX(0.D0,0.D0)

         CALL ZGEMM('n','n',N,P,M,ALP,A,LDA,B,LDB,BET,C,LDC)


         ! WRITE(6,*) 'In complex', N,M,P
         ! DO I = 1,N
         ! DO J = 1,P
         ! X = CMPLX(0.D0,0.D0)
         ! DO K = 1,M
         ! X = X + A(I,K)*B(K,J)
         ! ENDDO
         ! C(I,J) = X
         ! ENDDO
         ! ENDDO
         
       END SUBROUTINE MMULT_C

!*********
       SUBROUTINE INITD_R(A,X)
         IMPLICIT NONE
         REAL (KIND=8), DIMENSION(:,:) :: A
         REAL (KIND=8) X
         INTEGER I,J, N, M

         N = SIZE(A,1)
         M = SIZE(A,2)

         ! WRITE(6,*) 'In Init1 real', N,M
         DO I = 1,N
         DO J = 1,M
            A(I,J) = 0.D0
         ENDDO
         ENDDO
         DO I = 1,N
            A(I,I) = X
         ENDDO
       END SUBROUTINE INITD_R

       SUBROUTINE INITD_C(A,X)
         IMPLICIT NONE
         COMPLEX (KIND=8), DIMENSION(:,:) :: A
         COMPLEX (KIND=8) X
         INTEGER I,J, N, M

         N = SIZE(A,1)
         M = SIZE(A,2)

! WRITE(6,*) 'In Init1 complex', N,M
         DO I = 1,N
         DO J = 1,M
            A(I,J) = CMPLX(0.D0,0.D0)
         ENDDO
         ENDDO
         DO I = 1,N
            A(I,I) = X
         ENDDO
       END SUBROUTINE INITD_C


!*************
       SUBROUTINE INV_R0(A,AINV,DET)
         IMPLICIT NONE
         REAL (KIND=8), DIMENSION(:,:) :: A,AINV
         REAL (KIND=8) :: DET
         INTEGER I,J, N, M

! Working space.
         REAL (KIND=8) :: DET1(2)
         REAL (KIND=8), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT
         INTEGER INFO, JOB, LDA

         LDA = SIZE(A,1)
! Working space.
         ALLOCATE ( IPVT(LDA) )
         ALLOCATE ( WORK(LDA) )


         DO I = 1,LDA
            DO J = 1,LDA
               AINV(J,I) = A(J,I)
            ENDDO
         ENDDO

! Linpack routines.

         CALL DGEFA(AINV,LDA,LDA,IPVT,INFO)
         JOB = 11
         CALL DGEDI(AINV,LDA,LDA,IPVT,DET1,WORK,JOB)

         !Write(6,*) 'In Inv_R0', DET1
         DET = DET1(1) * 10.D0**DET1(2)


         DEALLOCATE (IPVT)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_R0


!*************
       SUBROUTINE INV_R_Variable(A,AINV,DET,Ndim)
         IMPLICIT NONE
         REAL (KIND=8), DIMENSION(:,:) :: A,AINV
         REAL (KIND=8) :: DET
         INTEGER I,J, N, M, Ndim

! Working space.
         REAL (KIND=8) :: DET1(2)
         REAL (KIND=8), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT
         INTEGER INFO, JOB, LDA

         LDA =  SIZE(A,1)
! Working space.
         ALLOCATE ( IPVT(Ndim) )
         ALLOCATE ( WORK(Ndim) )


         DO I = 1,LDA
            DO J = 1,LDA
               AINV(J,I) = A(J,I)
            ENDDO
         ENDDO

! Linpack routines.

         CALL DGEFA(AINV,LDA,Ndim,IPVT,INFO)
         JOB = 11
         CALL DGEDI(AINV,LDA,Ndim,IPVT,DET1,WORK,JOB)

         DET = DET1(1) * 10.D0**DET1(2)

         DEALLOCATE (IPVT)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_R_VARIABLE

!*************
       SUBROUTINE INV_R_Variable_1(A,AINV,DET,Ndim)
         IMPLICIT NONE
         REAL (KIND=8), DIMENSION(:,:) :: A,AINV
         REAL (KIND=8) :: DET(2)
         INTEGER I,J, N, M, Ndim

! Working space.
         REAL (KIND=8) :: DET1(2)
         REAL (KIND=8), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT
         INTEGER INFO, JOB, LDA

         LDA =  SIZE(A,1)
! Working space.
         ALLOCATE ( IPVT(Ndim) )
         ALLOCATE ( WORK(Ndim) )


         DO I = 1,LDA
            DO J = 1,LDA
               AINV(J,I) = A(J,I)
            ENDDO
         ENDDO

! Linpack routines.

         CALL DGEFA(AINV,LDA,Ndim,IPVT,INFO)
         JOB = 11
         CALL DGEDI(AINV,LDA,Ndim,IPVT,DET1,WORK,JOB)

         ! Determinant  = DET1(1) * 10.D0**DET1(2)
         DET(1) = DET1(1) 
         DET(2) = DET1(2)

         DEALLOCATE (IPVT)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_R_VARIABLE_1


!*************
       SUBROUTINE INV_R1(A,AINV,DET1)
         IMPLICIT NONE
         REAL (KIND=8), DIMENSION(:,:) :: A,AINV
         REAL (KIND=8) :: DET1(2)
         INTEGER I,J, N, M

! Working space.
         REAL (KIND=8), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT
         INTEGER INFO, JOB, LDA

         LDA = SIZE(A,1)
! Working space.
         ALLOCATE ( IPVT(LDA) )
         ALLOCATE ( WORK(LDA) )


         DO I = 1,LDA
            DO J = 1,LDA
               AINV(J,I) = A(J,I)
            ENDDO
         ENDDO

! Linpack routines.

         CALL DGEFA(AINV,LDA,LDA,IPVT,INFO)
         JOB = 11
         CALL DGEDI(AINV,LDA,LDA,IPVT,DET1,WORK,JOB)







         DEALLOCATE (IPVT)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_R1

!*************
       SUBROUTINE INV_R2(A,AINV)
         IMPLICIT NONE
         REAL (KIND=8), DIMENSION(:,:) :: A,AINV

         INTEGER I,J, N, M

! Uses Lapack routines.

! Working space.
         REAL (KIND=8), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
         INTEGER INFO, JOB, LDA, LWORK

         LDA = SIZE(A,1)

         !Write(6,*) 'Inv_r2:', LDA
         ALLOCATE ( IPIV(LDA) )
         LWORK = LDA
         ALLOCATE ( WORK(LWORK) )
         WORK = 0.0
         IPIV = 0
         DO I = 1,LDA
            DO J = 1,LDA
               AINV(J,I) = A(J,I)
            ENDDO
         ENDDO
         INFO = 0


         CALL DGETRF( LDA, LDA, AINV, LDA, IPIV, INFO )
         CALL DGETRI(LDA, AINV, LDA, IPIV, WORK, LWORK, INFO)







! Compute the determinant here if needed.
! detz = dcmplx(1.d0,0.d0)
! do n = 1,ne
! detz = detz * AINV(n,n)
! enddo ! Check. This may be wrong.


         DEALLOCATE (IPIV)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_R2
!*************

       SUBROUTINE INV_C(A,AINV,DET)
         IMPLICIT NONE
         COMPLEX (KIND=8), DIMENSION(:,:) :: A,AINV
         COMPLEX (KIND=8) :: DET
         INTEGER I,J, N, M

! Working space.
         COMPLEX (KIND=8) :: DET1(2)
         COMPLEX (KIND=8), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT
         INTEGER INFO, JOB, LDA

         LDA = SIZE(A,1)
! Working space.
         ALLOCATE ( IPVT(LDA) )
         ALLOCATE ( WORK(LDA) )


         DO I = 1,LDA
            DO J = 1,LDA
               AINV(J,I) = A(J,I)
            ENDDO
         ENDDO

! Linpack routines.

         CALL ZGEFA(AINV,LDA,LDA,IPVT,INFO)
         JOB = 11
         CALL ZGEDI(AINV,LDA,LDA,IPVT,DET1,WORK,JOB)







         DET = DET1(1)*10.D0**DET1(2)

         DEALLOCATE (IPVT)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_C

!========================================================================
       SUBROUTINE INV_C_Variable(A,AINV,DET,Ndim)
         IMPLICIT NONE
         COMPLEX (KIND=8), DIMENSION(:,:) :: A,AINV
         COMPLEX (KIND=8) :: DET
         INTEGER I,J, N, M,Ndim

! Working space.
         COMPLEX (KIND=8) :: DET1(2)
         COMPLEX (KIND=8), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT
         INTEGER INFO, JOB, LDA

         LDA = SIZE(A,1)
! Working space.
         ALLOCATE ( IPVT(Ndim) )
         ALLOCATE ( WORK(Ndim) )


         DO I = 1,LDA
            DO J = 1,LDA
               AINV(J,I) = A(J,I)
            ENDDO
         ENDDO

! Linpack routines.

         CALL ZGEFA(AINV,LDA,Ndim,IPVT,INFO)
         JOB = 11
         CALL ZGEDI(AINV,LDA,Ndim,IPVT,DET1,WORK,JOB)


         DET = DET1(1)*10.D0**DET1(2)

         DEALLOCATE (IPVT)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_C_VARIABLE

!========================================================================
       SUBROUTINE INV_C1(A,AINV,DET1)
         IMPLICIT NONE
         COMPLEX (KIND=8), DIMENSION(:,:) :: A,AINV
         COMPLEX (KIND=8) :: DET1(2)
         INTEGER I,J, N, M

! Working space.
         COMPLEX (KIND=8), DIMENSION(:), ALLOCATABLE :: WORK
         INTEGER, DIMENSION(:), ALLOCATABLE :: IPVT
         INTEGER INFO, JOB, LDA

         LDA = SIZE(A,1)
! Working space.
         ALLOCATE ( IPVT(LDA) )
         ALLOCATE ( WORK(LDA) )


         DO I = 1,LDA
            DO J = 1,LDA
               AINV(J,I) = A(J,I)
            ENDDO
         ENDDO

! Linpack routines.

         CALL ZGEFA(AINV,LDA,LDA,IPVT,INFO)
         JOB = 11
         CALL ZGEDI(AINV,LDA,LDA,IPVT,DET1,WORK,JOB)


         DEALLOCATE (IPVT)
         DEALLOCATE (WORK)
       END SUBROUTINE INV_C1
!*****



       SUBROUTINE COMPARE_C(A,B,XMAX,XMEAN)
         IMPLICIT NONE
         COMPLEX (KIND=8), DIMENSION(:,:) :: A,B
         REAL (KIND=8) :: XMAX, XMEAN
         INTEGER I,J, N, M

         REAL (KIND=8) :: DIFF

         N = SIZE(A,1)
         M = SIZE(A,2)

         XMAX = 0.D0
         XMEAN = 0.D0
         DO I = 1,N
            DO J = 1,M
               DIFF = SQRT( (A(I,J) - B(I,J))*CONJG(A(I,J)-B(I,J)))
               IF (DIFF.GT.XMAX) XMAX = DIFF
               XMEAN = XMEAN + DIFF
            ENDDO
         ENDDO
         XMEAN = XMEAN/DBLE(N*M)
       END SUBROUTINE COMPARE_C

       SUBROUTINE COMPARE_R(A,B,XMAX,XMEAN)
         IMPLICIT NONE
         REAL (KIND=8) , INTENT(IN), DIMENSION(:,:) :: A,B
         REAL (KIND=8) , INTENT(INOUT) :: XMAX, XMEAN
         INTEGER I,J, N, M

         REAL (KIND=8) :: DIFF

         N = SIZE(A,1)
         M = SIZE(A,2)

         XMAX = 0.D0
         XMEAN = 0.D0
         DO I = 1,N
            DO J = 1,M
               DIFF = ABS( ( B(I,J) - A(I,J) ) )
               IF (DIFF.GT.XMAX) XMAX = DIFF
               XMEAN = XMEAN + DIFF
            ENDDO
         ENDDO
         XMEAN = XMEAN/DBLE(N*M)
       END SUBROUTINE COMPARE_R

!*****************
       SUBROUTINE UDV1_R(A,U,D,V,NCON)
         IMPLICIT NONE
         REAL (KIND=8), INTENT(IN), DIMENSION(:,:) :: A
         REAL (KIND=8), INTENT(INOUT), DIMENSION(:,:) :: U,V
         REAL (KIND=8), INTENT(INOUT), DIMENSION(:) :: D
         INTEGER, INTENT(IN) :: NCON
         INTEGER I,J,K, N, M, ND1, ND2, NR, IMAX, IFAIL

!        The Det of V is not equal to unity. 
! Locals:
         INTEGER, DIMENSION(:), ALLOCATABLE :: IVPT, IVPTM1
         REAL (KIND=8), DIMENSION(:), ALLOCATABLE :: XNORM, VHELP,&
              & THETA, WORK
         REAL (KIND=8), DIMENSION(:,:), ALLOCATABLE :: TMP, V1,&
              & TEST, TEST1, TEST2
         REAL (KIND=8) :: XMAX, XMEAN, Z, DETV

         ND1 = SIZE(A,1)
         ND2 = SIZE(A,2)



! WRITE(6,*) 'Udv A: ',ND1,ND2
! WRITE(6,*) 'Udv V: ',size(V,1), size(V,2)
! You should now check corresponding sizes for U,V,D.
         IF (SIZE(U,1).NE.ND1 .OR. SIZE(U,2).NE.ND2) THEN
            WRITE(6,*) 'UDV dim mistake: U'
            STOP
         ENDIF
         IF (SIZE(D,1).NE.ND2 ) THEN
            WRITE(6,*) 'UDV dim mistake: D'
            STOP
         ENDIF
         IF (SIZE(V,1).NE.ND2 .OR. SIZE(V,2).NE.ND2) THEN
            WRITE(6,*) 'UDV dim mistake: V'
            STOP
         ENDIF

         ALLOCATE(XNORM (ND2))
         ALLOCATE(VHELP (ND2))
         ALLOCATE(IVPT (ND2))
         ALLOCATE(IVPTM1(ND2))
         ALLOCATE(WORK (ND2))
         ALLOCATE(THETA (ND2))

         ALLOCATE(TMP(ND1,ND2))
         ALLOCATE(V1 (ND2,ND2))

         V1 = 0.D0

         DO I = 1,ND2
            XNORM(I) = 0.D0
            DO NR = 1,ND1
               XNORM(I) = XNORM(I) + ABS(A(NR,I))
            ENDDO
         ENDDO
         DO I = 1,ND2
            VHELP(I) = XNORM(I)
         ENDDO

         DO I = 1,ND2
            XMAX = 0.D0
            DO J = 1,ND2
               IF (VHELP(J).GT.XMAX) IMAX = J
               IF (VHELP(J).GT.XMAX) XMAX = VHELP(J)
            ENDDO
            VHELP(IMAX) = -1.D0
            IVPTM1(IMAX)=I
            IVPT(I) = IMAX
         ENDDO

         DO I = 1,ND2
            Z = 1.D0/XNORM(IVPT(I))
            K = IVPT(I)
            DO NR = 1,ND1
               TMP(NR,I) = A(NR,K)*Z
            ENDDO
         ENDDO


         !You now want to UDV TMP. Nag routines.
         IFAIL = 0


         CALL F01QCF(ND1,ND2,TMP,ND1,THETA,IFAIL)
         

         !Scale V1 to a unit triangluar matrix.
         DO I = 1,ND2
            D(I) = ABS(TMP(I,I))
         ENDDO
         DO I = 1,ND2
            Z = 1.D0/D(I)
            DO J = I,ND2
               V1(I,J) = TMP(I,J)*Z
            ENDDO
         ENDDO
         
         
! Compute U
         IFAIL = 0

         CALL F01QEF('Separate', ND1,ND2, ND2, TMP,&
              & ND1, THETA, WORK, IFAIL)


         DO I = 1,ND1
            DO J = 1,ND2
               U(I,J) = TMP(I,J)
            ENDDO
         ENDDO


! Finish the pivotting.
         DO I = 1,ND2
            D(I) = D(I)*XNORM(IVPT(I))
         ENDDO
         DO I = 1,ND2-1
            Z = 1.D0/XNORM(IVPT(I))
            DO J = I+1,ND2
               V1(I,J) = V1(I,J)*XNORM(IVPT(J))*Z
            ENDDO
         ENDDO

         DO J = 1,ND2
            DO I = 1,ND2
               V(I,J) = V1(I,IVPTM1(J))
            ENDDO
         ENDDO

! Test accuracy.
         IF (NCON.EQ.1) THEN
            ALLOCATE (TEST(ND1,ND2))
            DO J = 1,ND2
               DO I = 1,ND1
                  Z = 0.D0
                  DO NR = 1,ND2
                     Z = Z + U(I,NR)*D(NR)*V(NR,J)
                  ENDDO
                  TEST(I,J) = Z
               ENDDO
            ENDDO
            XMAX = 0.0; XMEAN = 0.0
            CALL COMPARE(TEST,A,XMAX,XMEAN)
            WRITE(6,*) 'Accuracy: ',XMAX
            DEALLOCATE (TEST)

            ALLOCATE (TEST (ND2,ND1))
            ALLOCATE (TEST1 (ND2,ND2))
            ALLOCATE (TEST2 (ND2,ND2))
            ! Check orthogonality of U
            DO I = 1,ND1
               DO J = 1,ND2
                  TEST(J,I) = U(I,J)
               ENDDO
            ENDDO
            CALL MMULT(TEST1,TEST,U)
            CALL INITD(TEST2,1.D0)
            XMAX = 0.0; XMEAN = 0.0
            CALL COMPARE(TEST1,TEST2,XMAX,XMEAN)
            WRITE(6,*) 'UDV1 orth U: ',XMAX
            DEALLOCATE (TEST )
            DEALLOCATE (TEST1 )
            DEALLOCATE (TEST2 )
         ENDIF


         DEALLOCATE(XNORM )
         DEALLOCATE(VHELP )
         DEALLOCATE(IVPT )
         DEALLOCATE(IVPTM1)
         DEALLOCATE(WORK )
         DEALLOCATE(THETA )

         DEALLOCATE(TMP)
         DEALLOCATE(V1 )

      END SUBROUTINE UDV1_R

!***************
      SUBROUTINE UDV_C(A,U,D,V,NCON)
        !Uses Nag library.
        !#include "machine"

        IMPLICIT NONE
        COMPLEX (KIND=8), INTENT(IN), DIMENSION(:,:) :: A
        COMPLEX (KIND=8), INTENT(INOUT), DIMENSION(:,:) :: U,V
        COMPLEX (KIND=8), INTENT(INOUT), DIMENSION(:) :: D
        INTEGER, INTENT(IN) :: NCON
        INTEGER :: NE, LQ, IFAIL, I, J, NR

        !Local
        COMPLEX (KIND=8), DIMENSION(:,:), ALLOCATABLE :: TMP, TEST
        COMPLEX (KIND=8), DIMENSION(:), ALLOCATABLE :: THETA, WORK
        COMPLEX (KIND=8) :: Z
        REAL (KIND=8) :: DETV, XMDIFF, X

        LQ = SIZE(A,1)
        NE = SIZE(A,2)

        U = DCMPLX(0.D0,0.D0) ; V = DCMPLX(0.D0,0.D0); D = DCMPLX(0.D0,0.D0)
        ALLOCATE (TMP(LQ,NE), THETA(NE), WORK(NE))

        TMP = A

        !You now want to UDV TMP. Nag routines.
        IFAIL = 0

        CALL F01RCF(LQ,NE,TMP,LQ,THETA,IFAIL)


        DO I = 1,NE
           DO J = I,NE
              V(I,J) = TMP(I,J)
           ENDDO
        ENDDO
        DETV = 1.D0
        !V is an NE by NE upper triangular matrix with real diagonal elements.
        DO I = 1,NE
           DETV = DETV * DBLE( TMP(I,I) )
        ENDDO

        !Compute U

        CALL F01REF('Separate', LQ,NE, NE, TMP, &
             & LQ, THETA, WORK, IFAIL)





        DO J = 1,NE
           DO I = 1,LQ
              U(I,J) = TMP(I,J)
           ENDDO
        ENDDO

        IF (DBLE(DETV).LT.0.D0) THEN
           DO I = 1,LQ
              U(I,1) = -U(I,1)
           ENDDO
           DO I = 1,NE
              V(1,I) = -V(1,I)
           ENDDO
        ENDIF

        !Scale V1 to a unit triangluar matrix.
        DO I = 1,NE
           D(I) = CMPLX(ABS(DBLE(V(I,I))),0.D0)
        ENDDO
        DO I = 1,NE
           Z = DCMPLX(1.D0,0.D0)/D(I)
           DO J = I,NE
              V(I,J) = V(I,J)*Z
           ENDDO
        ENDDO

        !Test accuracy.
        IF (NCON.EQ.1) THEN
           ALLOCATE( TEST(LQ,NE) )
           DO J = 1,NE
              DO I = 1,LQ
                 Z = DCMPLX(0.D0,0.D0)
                 DO NR = 1,NE
                    Z = Z + U(I,NR)*D(NR)*V(NR,J)
                 ENDDO
                 TEST(I,J) = Z
              ENDDO
           ENDDO
           XMDIFF = 0.D0
           DO J = 1,LQ
              DO I = 1,NE
                 Z = (TEST(J,I)-A(J,I)) * CONJG(TEST(J,I)-A(J,I))
                 X = SQRT(DBLE(Z))
                 IF (X.GT.XMDIFF) XMDIFF = X
              ENDDO
           ENDDO
           WRITE(6,*) 'Accuracy, ortho: ',XMDIFF
           DEALLOCATE( TEST )
        ENDIF

        DEALLOCATE (TMP, THETA, WORK)

        RETURN
      END SUBROUTINE UDV_C

!***************
      SUBROUTINE QR_C(A,U,V,NCON)
        !Uses Nag library.
        !#include "machine"

        IMPLICIT NONE
        COMPLEX (KIND=8), INTENT(IN), DIMENSION(:,:) :: A
        COMPLEX (KIND=8), INTENT(INOUT), DIMENSION(:,:) :: U,V
        INTEGER, INTENT(IN) :: NCON
        INTEGER :: NE, LQ, IFAIL, I, J, NR

        !Local
        COMPLEX (KIND=8), DIMENSION(:,:), ALLOCATABLE :: TMP, TEST
        COMPLEX (KIND=8), DIMENSION(:), ALLOCATABLE :: THETA, WORK
        COMPLEX (KIND=8) :: Z
        REAL (KIND=8) :: DETV, XMDIFF, X

        LQ = SIZE(A,1)
        NE = SIZE(A,2)

        U = DCMPLX(0.D0,0.D0) ; V = DCMPLX(0.D0,0.D0)
        ALLOCATE (TMP(LQ,NE), THETA(NE), WORK(NE))

        TMP = A

        !You now want to UDV TMP. Nag routines.
        IFAIL = 0

        CALL F01RCF(LQ,NE,TMP,LQ,THETA,IFAIL)


        DO I = 1,NE
           DO J = I,NE
              V(I,J) = TMP(I,J)
           ENDDO
        ENDDO
        DETV = 1.D0
        !V is an NE by NE upper triangular matrix with real diagonal elements.
        DO I = 1,NE
           DETV = DETV * DBLE( TMP(I,I) )
        ENDDO

        !Compute U

        CALL F01REF('Separate', LQ,NE, NE, TMP, &
             & LQ, THETA, WORK, IFAIL)

        DO J = 1,NE
           DO I = 1,LQ
              U(I,J) = TMP(I,J)
           ENDDO
        ENDDO

        IF (DBLE(DETV).LT.0.D0) THEN
           DO I = 1,LQ
              U(I,1) = -U(I,1)
           ENDDO
           DO I = 1,NE
              V(1,I) = -V(1,I)
           ENDDO
        ENDIF


        !Test accuracy.
        IF (NCON.EQ.1) THEN
           ALLOCATE( TEST(LQ,NE) )
           DO J = 1,NE
              DO I = 1,LQ
                 Z = DCMPLX(0.D0,0.D0)
                 DO NR = 1,NE
                    Z = Z + U(I,NR)*V(NR,J)
                 ENDDO
                 TEST(I,J) = Z
              ENDDO
           ENDDO
           XMDIFF = 0.D0
           DO J = 1,LQ
              DO I = 1,NE
                 Z = (TEST(J,I)-A(J,I)) * CONJG(TEST(J,I)-A(J,I))
                 X = SQRT(DBLE(Z))
                 IF (X.GT.XMDIFF) XMDIFF = X
              ENDDO
           ENDDO
           WRITE(6,*) 'Accuracy, QR: ',XMDIFF
           DEALLOCATE( TEST )
        ENDIF

        DEALLOCATE (TMP, THETA, WORK)

        RETURN
      END SUBROUTINE QR_C
!********************
      SUBROUTINE SVD_C(A,U,D,V,NCON)
        !Uses LaPack Routine 
        !#include "machine"

        IMPLICIT NONE
        COMPLEX (KIND=8), INTENT(IN), DIMENSION(:,:) :: A
        COMPLEX (KIND=8), INTENT(INOUT), DIMENSION(:,:) :: U,V
        COMPLEX (KIND=8), INTENT(INOUT), DIMENSION(:) :: D
        INTEGER, INTENT(IN) :: NCON

        !! Local
        REAL    (Kind=8), Allocatable :: RWORK(:), S(:)
        COMPLEX (Kind=8), Allocatable :: WORK(:), A1(:,:)
        CHARACTER (Len=1):: JOBU,JOBVT
        INTEGER          :: M,N, LDA, LDVT, LDU, LWORK, I, J, I1, INFO
        REAL    (Kind=8) :: X, Xmax
        COMPLEX (Kind=8) :: Z
        
        JOBU = "A"
        JOBVT= "A"
        M = SIZE(A,1)
        N = SIZE(A,2)
        Allocate (A1(M,N))
        Allocate (S(N))
        A1  = A
        LDA = M
        LDU = M
        LDVT = N
        if (M > N) then
           LWORK = 2*N + M
           I = 3*N
           IF ( 5*N -4  > I) I =  5*N -4
           ALLOCATE (RWORK(I))
        Else
           LWORK = 2*M + N
           I = 3*M
           IF ( 5*M -4  > I) I =  5*M -4
           ALLOCATE (RWORK(I))
        Endif
        Allocate (WORK(LWORK))
        

        CALL ZGESVD( JOBU, JOBVT, M, N, A1, LDA, S, U, LDU, V, LDVT,&
             &        WORK, LWORK, RWORK, INFO )

        DO I = 1,N
           D(I) = cmplx(S(I),0.d0,kind=8)
        ENDDO

        IF (NCON ==  1) THEN
           Write(6,*) JobU, JobVT
           Xmax = 0.d0
           DO I = 1,M
              DO I1 = 1,N
                 Z = cmplx(0.d0,0.d0,Kind=8)
                 DO J = 1,N
                    Z  =  Z + U(I,J) *D(J) *V(J,I1)
                 ENDDO
                 X = sqrt(Real((Z - A(I,I1))*Conjg(Z - A(I,I1))))
                 IF (X > Xmax ) Xmax = X
              ENDDO
           ENDDO
           WRITE(6,*) "Success (0), PRE ", INFO, Xmax
        ENDIF

        
        Deallocate (WORK,RWORK,A1,S)
        

      END SUBROUTINE SVD_C
!***************

      SUBROUTINE DIAG_R(A,U,W)
        IMPLICIT NONE
        REAL (KIND=8), INTENT(IN), DIMENSION(:,:) :: A
        REAL (KIND=8), INTENT(INOUT), DIMENSION(:,:) :: U
        REAL (KIND=8), INTENT(INOUT), DIMENSION(:) :: W


        INTEGER ND1,ND2,IERR
        REAL (KIND=8), DIMENSION(:), ALLOCATABLE :: WORK

         ND1 = SIZE(A,1)
         ND2 = SIZE(A,2)

         IF (ND1.NE.ND2) THEN
            WRITE(6,*) 'Error in matrix dimension DIAG_R'
            STOP
         ENDIF

         IERR = 0
         U=A
         W=0
         ! let's just give lapack enough memory
         ALLOCATE(WORK(3*ND1))
         CALL DSYEV('V', 'U', ND1, U, ND1, W, WORK, 3*ND1, IERR)
         DEALLOCATE(WORK)

      END SUBROUTINE DIAG_R
!*********

      SUBROUTINE DIAG_I(A,U,W)
        ! Uses Lapack
        IMPLICIT NONE
        COMPLEX (KIND=8), INTENT(IN)   , DIMENSION(:,:) :: A
        COMPLEX (KIND=8), INTENT(INOUT), DIMENSION(:,:) :: U
        REAL    (KIND=8), INTENT(INOUT), DIMENSION(:)   :: W
        
        CHARACTER (len=1) :: UPLO, JOBZ
        INTEGER :: N, LWORK, INFO
        COMPLEX (KIND=8), allocatable :: WORK (:)
        REAL    (KIND=8), allocatable :: RWORK(:)
        Logical :: Test
        Integer :: I,J,m
        Complex (Kind=8) :: Z 
        Real (Kind=8) :: X, XMAX

        JOBZ = "V"
        UPLO = "U"
        N = size(A,1) 
        U = A
        LWORK = 2*N -1
        Allocate ( WORK(LWORK) )
        Allocate (  RWORK(3*N-2))
	
	!Write(6,*) 'In Diag'

        Call ZHEEV (JOBZ, UPLO, N, U, N, W, WORK, LWORK, RWORK, INFO)

        Deallocate (WORK, RWORK)
        
        Test = .false.
        If (Test) then
           XMAX = 0.d0
           DO I = 1,N
              DO J = 1,N
                 Z = cmplx(0.d0,0.d0,kind=8)
                 DO m = 1,N
                    Z =  Z + U(I,m)*cmplx(W(m),0.d0, Kind=8)*Conjg(U(J,m))
                 ENDDO
                 Z = Z - A(I,J)
                 X = sqrt( Z*Conjg(Z) )
                 If (X > XMAX ) XMAX = X
              ENDDO
           ENDDO
           write(6,*) ' Test Diag_I: ', XMAX
        endif
        
      End SUBROUTINE DIAG_I

      SUBROUTINE SECONDS(X)
        IMPLICIT NONE
        REAL (KIND=8), INTENT(INOUT) :: X

        !DATE_AND_TIME(date, time, zone, values)
        !date_and_time([date][,time][,zone][,values])
        !Subroutine. Die Parameter haben das Attribut intent(out), geben also Werte zur�ck.

        ! date: skalare, normale Zeichenvariable von wenigstens 8 Zeichen. Die linken 8 Zeichen bekommen einen Wert der Form JJJJMMTT . JJJJ Jahr, MM Monat, TT Tag im Monat.
        !time: skalare, normale Zeichenvariable von wenigstens 10 Zeichen. Die linken 10 Zeichen bekommen einen Wert der Form hhmmss.sss , wobei hh die Stunde des Tages ist, mm die Minute innerhalb der Stunde, und ss.sss die Sekunde mit Bruchteilen.
        ! zone: skalare, normale Zeichenvariable von wenigstens 5 Zeichen. Die linken 5 Zeichen bekommen einen Wert der Form hhmm . hh Stunden, mm Minuten Zeitdifferenz gegen�ber der UTC-Weltzeit.
        !values: Eindimensionales Integer-Feld. L�nge wenigstens 8. 1: Jahr, z.B. 1993. 2: Monat. 3: Monatstag. 4: Zeitdifferenz zur Weltzeit in Minuten. 5: Stunde des Tages. 6: Minute innerhalb der Stunde. 7: Sekunden 8. Millisekunden.

        !character(len=10) :: d,t
        integer,dimension(8) :: V
        !d = ""
        !call date_and_time(date=d,time=t)
        call date_and_time(values=V)

        X = DBLE(V(5)*3600 + V(6)*60 + V(7))

      END SUBROUTINE SECONDS

!====================================================
      Complex (Kind=8) Function DET_C(Mat,N)

        Implicit none
        
        ! Arguments
        Integer, intent(in) :: N 
        Complex(kind=8), intent(inout) :: mat(N,N)
        
        integer :: i, info
        integer :: ipiv(N)

        integer :: sgn

        ipiv = 0

        !Lapack LU decomposition
        call zgetrf(N, N, mat, N, ipiv, info)
        
        det_C = cmplx(1.d0,0.d0) 
        do i = 1, N
           det_C = det_C*mat(i, i)
        enddo
        
        sgn =  1
        do i = 1, N
           if(ipiv(i) /= i)  sgn = -sgn
        enddo
        if (sgn == -1 ) det_C = - det_C 

      end function DET_C
        

    END MODULE MyMats
