subroutine ZSLMM(side, op, N, M, A, P, Mat)

        IMPLICIT NONE
        CHARACTER (1)            , INTENT(IN) :: side, op
        INTEGER                  , INTENT(IN) :: N, M
        COMPLEX (KIND=KIND(0.D0)), INTENT(IN)   , DIMENSION(N,N) :: A
        COMPLEX (KIND=KIND(0.D0)), INTENT(INOUT), DIMENSION(M,M) :: Mat
        INTEGER                  , INTENT(IN)   , DIMENSION(N)   :: P
        
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: WORK, WORK2
        Complex (Kind = Kind(0.D0)) :: alpha, beta, Z(4)
        INTEGER :: I,L,IDX, NUMBLOCKS
        INTEGER, DIMENSION(:), ALLOCATABLE :: IDXLIST, DIMLIST
        LOGICAL :: COMPACT, LEFT
        
        
        IF ( side == 'L' .or. side == 'l' ) THEN
          ! multiply op(A) from the left [ Mat = op(A)*Mat ]
          LEFT = .TRUE.
        ELSEIF ( side == 'R' .or. side == 'r' ) THEN
          LEFT = .FALSE.
        ELSE
          write(*,*) 'Illegal argument for side: It is not one of [R,r,L,l] !'
          call EXIT(1)
        ENDIF
        
        alpha = 1.D0
        beta = 0.D0
        
        !identify possible block structure
        COMPACT = .TRUE.
        L = 1
        IDX = 1
        ALLOCATE(IDXLIST(N),DIMLIST(N))
        NUMBLOCKS=0
        DO I=1,N-1
          IF ( P(I)+1 .ne. P(I+1) ) THEN
            COMPACT = .FALSE.
            NUMBLOCKS=NUMBLOCKS+1
            IDXLIST(NUMBLOCKS)=IDX
            DIMLIST(NUMBLOCKS)=L
            IDX=IDX+L
            L=1
          ELSE
            L=L+1
          ENDIF
        ENDDO
        IF(IDX<N+1) THEN
          ! last block if need be
            NUMBLOCKS=NUMBLOCKS+1
            IDXLIST(NUMBLOCKS)=IDX
            DIMLIST(NUMBLOCKS)=L
        ENDIF
        
        IF (LEFT) THEN
          ! multiply op(A) from the left [ Mat = op(A)*Mat ]
          SELECT CASE(N)
          CASE (1)
            ! Here only one row is rescaled
            CALL ZSCAL(M,A(1,1),Mat(P(1),1),M)
          CASE (2)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Mat(P(1),I)=A(1,1)*Z(1)+A(1,2)*Z(2)
              Mat(P(2),I)=A(2,1)*Z(1)+A(2,2)*Z(2)
            ENDDO
          CASE (3)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Mat(P(1),I)=A(1,1)*Z(1)+A(1,2)*Z(2)+A(1,3)*Z(3)
              Mat(P(2),I)=A(2,1)*Z(1)+A(2,2)*Z(2)+A(2,3)*Z(3)
              Mat(P(3),I)=A(3,1)*Z(1)+A(3,2)*Z(2)+A(3,3)*Z(3)
            ENDDO
          CASE (4)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Mat(P(1),I)=A(1,1)*Z(1)+A(1,2)*Z(2)+A(1,3)*Z(3)+A(1,4)*Z(4)
              Mat(P(2),I)=A(2,1)*Z(1)+A(2,2)*Z(2)+A(2,3)*Z(3)+A(2,4)*Z(4)
              Mat(P(3),I)=A(3,1)*Z(1)+A(3,2)*Z(2)+A(3,3)*Z(3)+A(3,4)*Z(4)
              Mat(P(4),I)=A(4,1)*Z(1)+A(4,2)*Z(2)+A(4,3)*Z(3)+A(4,4)*Z(4)
            ENDDO
          CASE DEFAULT
            ! allocate memory and copy blocks of Mat to work
            ALLOCATE(WORK(N,M))
            DO I=1,NUMBLOCKS
              CALL ZLACPY('A', DIMLIST(I), M, Mat(P(IDXLIST(I)),1), M, WORK(IDXLIST(I),1), N)
            ENDDO
            
            ! Perform Mat multiplication
            IF(COMPACT) THEN
              !write result directly into mat
              CALL ZGEMM(op,'N', N, M, N, alpha, A(1, 1), N, WORK(1, 1), N, beta, Mat(P(1), 1), M)
            ELSE
              !additional space for result
              ALLOCATE(WORK2(N,M))
              CALL ZGEMM(op,'N', N, M, N, alpha, A(1, 1), N, WORK(1, 1), N, beta, WORK2(1, 1), N)
              !distribute result back into mat using blocks
              DO I=1,NUMBLOCKS
                CALL ZLACPY('A', DIMLIST(I), M, WORK2(IDXLIST(I),1), N, Mat(P(IDXLIST(I)),1), M)
              ENDDO
              !free result memory
              DEALLOCATE(WORK2)
            ENDIF
            !free memory of first mat copy
            DEALLOCATE(WORK)
          END SELECT
        ELSE
          ! multiply op(A) from the right [ Mat = Mat*op(A) ]
          SELECT CASE(N)
          CASE (1)
            ! Here only one column is rescaled
            CALL ZSCAL(M,A(1,1),Mat(1,P(1)),1)
          CASE (2)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Mat(I,P(1))=A(1,1)*Z(1)+A(2,1)*Z(2)
              Mat(I,P(2))=A(1,2)*Z(1)+A(2,2)*Z(2)
            ENDDO
          CASE (3)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Mat(I,P(1))=A(1,1)*Z(1)+A(2,1)*Z(2)+A(3,1)*Z(3)
              Mat(I,P(2))=A(1,2)*Z(1)+A(2,2)*Z(2)+A(3,2)*Z(3)
              Mat(I,P(3))=A(1,3)*Z(1)+A(2,3)*Z(2)+A(3,3)*Z(3)
            ENDDO
          CASE (4)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Mat(I,P(1))=A(1,1)*Z(1)+A(2,1)*Z(2)+A(3,1)*Z(3)+A(4,1)*Z(4)
              Mat(I,P(2))=A(1,2)*Z(1)+A(2,2)*Z(2)+A(3,2)*Z(3)+A(4,2)*Z(4)
              Mat(I,P(3))=A(1,3)*Z(1)+A(2,3)*Z(2)+A(3,3)*Z(3)+A(4,3)*Z(4)
              Mat(I,P(4))=A(1,4)*Z(1)+A(2,4)*Z(2)+A(3,4)*Z(3)+A(4,4)*Z(4)
            ENDDO
          CASE DEFAULT
            ! allocate memory and copy blocks of Mat to work
            ALLOCATE(WORK(M,N))
            DO I=1,NUMBLOCKS
              CALL ZLACPY('A', M, DIMLIST(I), Mat(1,P(IDXLIST(I))), M, WORK(1,IDXLIST(I)), M)
            ENDDO
            
            ! Perform Mat multiplication
            IF(COMPACT) THEN
              !write result directly into mat
              CALL ZGEMM('N',op, M, N, N, alpha, WORK(1, 1), M, A(1, 1), N, beta, Mat(1, P(1)), M)
            ELSE
              !additional space for result
              ALLOCATE(WORK2(M,N))
              CALL ZGEMM('N',op, M, N, N, alpha, WORK(1, 1), M, A(1, 1), N, beta, WORK2(1, 1), M)
              !distribute result back into mat using blocks
              DO I=1,NUMBLOCKS
                CALL ZLACPY('A', M, DIMLIST(I), WORK2(1,IDXLIST(I)), M, Mat(1,P(IDXLIST(I))), M)
              ENDDO
              !free result memory
              DEALLOCATE(WORK2)
            ENDIF 
            !free memory of first mat copy
            DEALLOCATE(WORK)
          END SELECT
        ENDIF
        DEALLOCATE(IDXLIST,DIMLIST)

end subroutine ZSLMM

PROGRAM Main

        IMPLICIT NONE
        INTEGER                   :: N, M, I, J, K, O
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: A
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: Mat, MatL, MatR
        INTEGER                  , DIMENSION(:)  , ALLOCATABLE :: P
        COMPLEX (KIND=KIND(0.D0)) :: Z
        
        M=400
        
        ALLOCATE(Mat(M,M),MatL(M,M),MatR(M,M))
        
        Do N=1,6
        
        ALLOCATE(A(N,N),P(N))
        
        Do O=7,M
        P(1)=O
        if (N>1) P(2)=2
        if (N>2) P(3)=1
        if (N>3) P(4)=3
        if (N>4) P(5)=5
        if (N>5) P(6)=6
        
        DO I=1,M
        DO J=1,M
          Mat(I,J) = CMPLX(DBLE(i+j),DBLE(i-j),kind(0.d0))
        ENDDO
        ENDDO
        
        DO I=1,N
        DO J=1,N
          A(I,J) = CMPLX(DBLE(200*i+10*j),DBLE(4*i+3*j),kind(0.d0))
        ENDDO
        ENDDO
        
        MatL(:,:)=Mat(:,:)
        DO I=1,N
        DO J=1,M
        Z = CMPLX(0.d0,0.d0,kind(0.d0))
        DO K=1,N
          Z = Z + A(I,K)*Mat(P(K),J) 
        ENDDO
        MatL(P(I),J) = Z
        ENDDO
        ENDDO
        
        MatR(:,:)=MatL(:,:)
        DO I=1,M
        DO J=1,N
        Z = CMPLX(0.d0,0.d0,kind(0.d0))
        DO K=1,N
          Z = Z + MatL(I,P(K)) * A(k,j)
        ENDDO
        MatR(I,P(J)) = Z
        ENDDO
        ENDDO
        
        
        CALL ZSLMM('L','N', N, M, A, P, Mat)
        IF(abs(sum(Mat-MatL))>1E-14)  write(*,*) "ERROR"
        
        CALL ZSLMM('R','N', N, M, A, P, Mat)
        IF(abs(sum(Mat-MatR))>1E-14)  write(*,*) "ERROR"
        
        enddo
        
        Deallocate(A,P)
        enddo
        Deallocate(Mat,MatL,MatR)
        write(*,*) "SUCCESS"
END PROGRAM Main

!     select case (opn)
!     case (1)
!         DO I = 1, Ndim
!             Mat(P(1), I) = V(1, I)
!         enddo
!     case (2)
!         DO I = 1, Ndim
!             Mat(P(1), I) = U(1, 1) * V(1, I) - conjg(U(2, 1)) * V(2, I)
!             Mat(P(2), I) = U(2, 1) * V(1, I) + conjg(U(1, 1)) * V(2, I)
!         enddo
!     case default
!         Allocate(tmp(opn, Ndim))
!         CALL ZGEMM('N','N', opn, Ndim, opn, alpha, U(1, 1), opn, V(1, 1), opn, beta, tmp(1, 1), opn)
!         Mat((P), :) = tmp
!         Deallocate(tmp)
!     end select
