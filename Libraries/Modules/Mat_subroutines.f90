subroutine ZSLGEMM(side, op, N, M, A, P, Mat)

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
        LOGICAL :: COMPACT, ALLOC
        
        alpha = 1.D0
        beta = 0.D0
        
        !identify possible block structure
        !only used in default case for n>4
        IF(N > 4) THEN
          COMPACT = .TRUE.
          L = 1
          IDX = 1
          ALLOCATE(IDXLIST(N),DIMLIST(N))
          ALLOC = .TRUE.
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
        ELSEIF(N>1) THEN
          ALLOCATE(WORK(N,N))
          IF( op == 'N' .or. op == 'n') THEN
            CALL ZLACPY('A', N, N, A(1,1), N, WORK(1,1), N)
          ELSEIF( op == 'T' .or. op == 't' ) then
            DO I=1,N
            DO L=1,N
              WORK(I,L)=A(L,I)
            ENDDO
            ENDDO
          ELSEIF( op == 'C' .or. op == 'c' ) then
            DO I=1,N
            DO L=1,N
              WORK(I,L)=conjg(A(L,I))
            ENDDO
            ENDDO
          ENDIF
        ENDIF
        
        IF ( side == 'L' .or. side == 'l' ) THEN
          ! multiply op(A) from the left  [ Mat = op(A)*Mat ]
          
          SELECT CASE(N)
          CASE (1)
            ! Here only one row is rescaled
            IF(op == 'N' .or. op=='n' .or. op=='T' .or. op=='t') then
              CALL ZSCAL(M,A(1,1),Mat(P(1),1),M)
            else
              CALL ZSCAL(M,conjg(A(1,1)),Mat(P(1),1),M)
            endif
          CASE (2)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)
            ENDDO
            DEALLOCATE(WORK)
          CASE (3)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)
            ENDDO
            DEALLOCATE(WORK)
          CASE (4)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)+WORK(1,4)*Z(4)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)+WORK(2,4)*Z(4)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)+WORK(3,4)*Z(4)
              Mat(P(4),I)=WORK(4,1)*Z(1)+WORK(4,2)*Z(2)+WORK(4,3)*Z(3)+WORK(4,4)*Z(4)
            ENDDO
            DEALLOCATE(WORK)
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
            DEALLOCATE(WORK,IDXLIST,DIMLIST)
          END SELECT
          
        ELSEIF ( side == 'R' .or. side == 'r' ) THEN
          ! multiply op(A) from the right [ Mat = Mat*op(A) ]
          
          SELECT CASE(N)
          CASE (1)
            ! Here only one column is rescaled
            IF(op == 'N' .or. op=='n' .or. op=='T' .or. op=='t') then
              CALL ZSCAL(M,A(1,1),Mat(1,P(1)),1)
            ELSE
              CALL ZSCAL(M,conjg(A(1,1)),Mat(1,P(1)),1)
            ENDIF
          CASE (2)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)
            ENDDO
            DEALLOCATE(WORK)
          CASE (3)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)
            ENDDO
            DEALLOCATE(WORK)
          CASE (4)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)+WORK(4,1)*Z(4)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)+WORK(4,2)*Z(4)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)+WORK(4,3)*Z(4)
              Mat(I,P(4))=WORK(1,4)*Z(1)+WORK(2,4)*Z(2)+WORK(3,4)*Z(3)+WORK(4,4)*Z(4)
            ENDDO
            DEALLOCATE(WORK)
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
            DEALLOCATE(WORK,IDXLIST,DIMLIST)
          END SELECT
          
        ELSE
          write(*,*) 'Illegal argument for side: It is not one of [R,r,L,l] !'
          call EXIT(1)
        ENDIF

end subroutine ZSLGEMM

subroutine ZSLHEMM(side, uplo, N, M, A, P, Mat)

        IMPLICIT NONE
        CHARACTER (1)            , INTENT(IN) :: side, uplo
        INTEGER                  , INTENT(IN) :: N, M
        COMPLEX (KIND=KIND(0.D0)), INTENT(IN)   , DIMENSION(N,N) :: A
        COMPLEX (KIND=KIND(0.D0)), INTENT(INOUT), DIMENSION(M,M) :: Mat
        INTEGER                  , INTENT(IN)   , DIMENSION(N)   :: P
        
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: WORK, WORK2
        Complex (Kind = Kind(0.D0)) :: alpha, beta, Z(4)
        INTEGER :: I,L,IDX, NUMBLOCKS
        INTEGER, DIMENSION(:), ALLOCATABLE :: IDXLIST, DIMLIST
        LOGICAL :: COMPACT, ALLOC
        
        alpha = 1.D0
        beta = 0.D0
        
        !identify possible block structure
        !only used in default case for n>4
        IF(N > 4) THEN
          COMPACT = .TRUE.
          L = 1
          IDX = 1
          ALLOCATE(IDXLIST(N),DIMLIST(N))
          ALLOC = .TRUE.
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
        ELSEIF(N>1) THEN
          ALLOCATE(WORK(N,N))
          CALL ZLACPY(uplo, N, N, A(1,1), N, WORK(1,1), N)
          ! Fill the rest of WORK (thereby ignoring uplo)
          IF(uplo=='U' .or. uplo=='u') THEN
            DO L=1,N
            DO I=1,L-1
              WORK(L,I)=conjg(WORK(I,L))
            ENDDO
            ENDDO
          ELSE
            DO L=1,N
            DO I=L+1,N
              WORK(L,I)=conjg(WORK(I,L))
            ENDDO
            ENDDO
          ENDIF
        ENDIF
        
        IF ( side == 'L' .or. side == 'l' ) THEN
          ! multiply op(A) from the left  [ Mat = op(A)*Mat ]
          
          SELECT CASE(N)
          CASE (1)
            ! Here only one row is rescaled
            ! uplo and transpositions as well as conjugatio has no effect for 1x1 herm. Matrices
            CALL ZSCAL(M,A(1,1),Mat(P(1),1),M)
          CASE (2)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)
            ENDDO
            DEALLOCATE(WORK)
          CASE (3)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)
            ENDDO
            DEALLOCATE(WORK)
          CASE (4)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)+WORK(1,4)*Z(4)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)+WORK(2,4)*Z(4)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)+WORK(3,4)*Z(4)
              Mat(P(4),I)=WORK(4,1)*Z(1)+WORK(4,2)*Z(2)+WORK(4,3)*Z(3)+WORK(4,4)*Z(4)
            ENDDO
            DEALLOCATE(WORK)
          CASE DEFAULT
            ! allocate memory and copy blocks of Mat to work
            ALLOCATE(WORK(N,M))
            DO I=1,NUMBLOCKS
              CALL ZLACPY('A', DIMLIST(I), M, Mat(P(IDXLIST(I)),1), M, WORK(IDXLIST(I),1), N)
            ENDDO
            
            ! Perform Mat multiplication
            IF(COMPACT) THEN
              !write result directly into mat
              CALL ZHEMM(side, uplo, N, M, alpha, A(1, 1), N, WORK(1, 1), N, beta, Mat(P(1), 1), M)
            ELSE
              !additional space for result
              ALLOCATE(WORK2(N,M))
              CALL ZHEMM(side, uplo, N, M, alpha, A(1, 1), N, WORK(1, 1), N, beta, WORK2(1, 1), N)
              !distribute result back into mat using blocks
              DO I=1,NUMBLOCKS
                CALL ZLACPY('A', DIMLIST(I), M, WORK2(IDXLIST(I),1), N, Mat(P(IDXLIST(I)),1), M)
              ENDDO
              !free result memory
              DEALLOCATE(WORK2)
            ENDIF
            !free memory of first mat copy
            DEALLOCATE(WORK,IDXLIST,DIMLIST)
          END SELECT
          
        ELSEIF ( side == 'R' .or. side == 'r' ) THEN
          ! multiply op(A) from the right [ Mat = Mat*op(A) ]
          
          SELECT CASE(N)
          CASE (1)
            ! Here only one column is rescaled
            ! uplo and transpositions as well as conjugatio has no effect for 1x1 herm. Matrices
            CALL ZSCAL(M,A(1,1),Mat(1,P(1)),1)
          CASE (2)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)
            ENDDO
            DEALLOCATE(WORK)
          CASE (3)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)
            ENDDO
            DEALLOCATE(WORK)
          CASE (4)
            ! perform inplace matmult
            DO I=1,M
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)+WORK(4,1)*Z(4)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)+WORK(4,2)*Z(4)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)+WORK(4,3)*Z(4)
              Mat(I,P(4))=WORK(1,4)*Z(1)+WORK(2,4)*Z(2)+WORK(3,4)*Z(3)+WORK(4,4)*Z(4)
            ENDDO
            DEALLOCATE(WORK)
          CASE DEFAULT
            ! allocate memory and copy blocks of Mat to work
            ALLOCATE(WORK(M,N))
            DO I=1,NUMBLOCKS
              CALL ZLACPY('A', M, DIMLIST(I), Mat(1,P(IDXLIST(I))), M, WORK(1,IDXLIST(I)), M)
            ENDDO
            
            ! Perform Mat multiplication
            IF(COMPACT) THEN
              !write result directly into mat
              CALL ZHEMM(side, uplo, M, N, alpha, A(1, 1), N, WORK(1, 1), M, beta, Mat(1, P(1)), M)
            ELSE
              !additional space for result
              ALLOCATE(WORK2(M,N))
              CALL ZHEMM(side, uplo, M, N, alpha, A(1, 1), N, WORK(1, 1), M, beta, WORK2(1, 1), M)
              !distribute result back into mat using blocks
              DO I=1,NUMBLOCKS
                CALL ZLACPY('A', M, DIMLIST(I), WORK2(1,IDXLIST(I)), M, Mat(1,P(IDXLIST(I))), M)
              ENDDO
              !free result memory
              DEALLOCATE(WORK2)
            ENDIF 
            !free memory of first mat copy
            DEALLOCATE(WORK,IDXLIST,DIMLIST)
          END SELECT
          
        ELSE
          write(*,*) 'Illegal argument for side: It is not one of [R,r,L,l] !'
          call EXIT(1)
        ENDIF

end subroutine ZSLHEMM

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
        
!         Do O=7,M
        P(1)=1
        if (N>1) P(2)=2
        if (N>2) P(3)=5
        if (N>3) P(4)=6
        if (N>4) P(5)=M-1
        if (N>5) P(6)=M
        
        DO I=1,N
        DO J=1,N
          A(I,J) = CMPLX(DBLE(i+j),DBLE(i-j),kind(0.d0))
        ENDDO
        ENDDO
        
        !test left-mult using upper
        DO I=1,M
        DO J=1,M
          Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
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
        
        CALL ZSLHEMM('L','u', N, M, A, P, Mat)
        IF(abs(sum(Mat-MatL))/M/M>1E-14) then
          write(*,*) "ERROR in ZSLHEMM left upper"
          write(*,*) abs(sum(Mat-MatL))/M/M, " is larger then 1E-14!"
          stop
        endif
        
        !test left-mult using lower
        DO I=1,M
        DO J=1,M
          Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
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
        
        CALL ZSLHEMM('L','l', N, M, A, P, Mat)
        IF(abs(sum(Mat-MatL))/M/M>1E-14) then
          write(*,*) "ERROR in ZSLHEMM left lower"
          write(*,*) abs(sum(Mat-MatL))/M/M, " is larger then 1E-14!"
          stop
        endif
        
        !test right-mult using upper
        DO I=1,M
        DO J=1,M
          Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
        ENDDO
        ENDDO
        MatR(:,:)=Mat(:,:)
        DO I=1,M
        DO J=1,N
        Z = CMPLX(0.d0,0.d0,kind(0.d0))
        DO K=1,N
          Z = Z + Mat(I,P(K)) * A(k,j)
        ENDDO
        MatR(I,P(J)) = Z
        ENDDO
        ENDDO
        
        CALL ZSLHEMM('R','u', N, M, A, P, Mat)
        IF(abs(sum(Mat-MatR))/M/M>1E-14) then
          write(*,*) "ERROR in ZSLHEMM right upper"
          write(*,*) abs(sum(Mat-MatR))/M/M, " is larger then 1E-14!"
          stop
        endif
        
        !test right-mult using lower
        DO I=1,M
        DO J=1,M
          Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
        ENDDO
        ENDDO
        MatR(:,:)=Mat(:,:)
        DO I=1,M
        DO J=1,N
        Z = CMPLX(0.d0,0.d0,kind(0.d0))
        DO K=1,N
          Z = Z + Mat(I,P(K)) * A(k,j)
        ENDDO
        MatR(I,P(J)) = Z
        ENDDO
        ENDDO
        
        CALL ZSLHEMM('R','l', N, M, A, P, Mat)
        IF(abs(sum(Mat-MatR))/M/M>1E-14) then
          write(*,*) "ERROR in ZSLHEMM right lower"
          write(*,*) abs(sum(Mat-MatR))/M/M, " is larger then 1E-14!"
          stop
        endif
        
        !test left-mult op='n'
        DO I=1,M
        DO J=1,M
          Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
        ENDDO
        ENDDO
        MatL(:,:)=Mat(:,:)
        DO I=1,N
        DO J=1,M
        Z = CMPLX(0.d0,0.d0,kind(0.d0))
        DO K=1,N
          Z = Z + A(i,k)*Mat(P(K),J) 
        ENDDO
        MatL(P(I),J) = Z
        ENDDO
        ENDDO
        
        CALL ZSLGEMM('L','n', N, M, A, P, Mat)
        IF(abs(sum(Mat-MatL))/M/M>1E-14) then
          write(*,*) "ERROR in ZSLGEMM left no-transpose"
          write(*,*) abs(sum(Mat-MatL))/M/M, " is larger then 1E-14!"
          stop
        endif
        
        !test left-mult op='t'
        DO I=1,M
        DO J=1,M
          Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
        ENDDO
        ENDDO
        MatL(:,:)=Mat(:,:)
        DO I=1,N
        DO J=1,M
        Z = CMPLX(0.d0,0.d0,kind(0.d0))
        DO K=1,N
          Z = Z + A(k,i)*Mat(P(K),J) 
        ENDDO
        MatL(P(I),J) = Z
        ENDDO
        ENDDO
        
        CALL ZSLGEMM('L','t', N, M, A, P, Mat)
        IF(abs(sum(Mat-MatL))/M/M>1E-14) then
          write(*,*) "ERROR in ZSLGEMM left transpose"
          write(*,*) abs(sum(Mat-MatL))/M/M, " is larger then 1E-14!"
          stop
        endif
        
        !test left-mult op='c'
        DO I=1,M
        DO J=1,M
          Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
        ENDDO
        ENDDO
        MatL(:,:)=Mat(:,:)
        DO I=1,N
        DO J=1,M
        Z = CMPLX(0.d0,0.d0,kind(0.d0))
        DO K=1,N
          Z = Z + conjg(A(k,i))*Mat(P(K),J) 
        ENDDO
        MatL(P(I),J) = Z
        ENDDO
        ENDDO
        
        CALL ZSLGEMM('L','c', N, M, A, P, Mat)
        IF(abs(sum(Mat-MatL))/M/M>1E-14) then
          write(*,*) "ERROR in ZSLGEMM left conjg-transpose"
          write(*,*) abs(sum(Mat-MatL))/M/M, " is larger then 1E-14!"
          stop
        endif
        
        !test right-mult op='n'
        DO I=1,M
        DO J=1,M
          Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
        ENDDO
        ENDDO
        MatR(:,:)=Mat(:,:)
        DO I=1,M
        DO J=1,N
        Z = CMPLX(0.d0,0.d0,kind(0.d0))
        DO K=1,N
          Z = Z + Mat(I,P(K)) * A(k,j)
        ENDDO
        MatR(I,P(J)) = Z
        ENDDO
        ENDDO
        
        CALL ZSLGEMM('R','n', N, M, A, P, Mat)
        IF(abs(sum(Mat-MatR))/M/M>1E-14) then
          write(*,*) "ERROR in ZSLGEMM right no-transpose"
          write(*,*) abs(sum(Mat-MatR))/M/M, " is larger then 1E-14!"
          stop
        endif
        
        !test right-mult op='t'
        DO I=1,M
        DO J=1,M
          Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
        ENDDO
        ENDDO
        MatR(:,:)=Mat(:,:)
        DO I=1,M
        DO J=1,N
        Z = CMPLX(0.d0,0.d0,kind(0.d0))
        DO K=1,N
          Z = Z + Mat(I,P(K)) * A(j,k)
        ENDDO
        MatR(I,P(J)) = Z
        ENDDO
        ENDDO
        
        CALL ZSLGEMM('R','T', N, M, A, P, Mat)
        IF(abs(sum(Mat-MatR))/M/M>1E-14) then
          write(*,*) "ERROR in ZSLGEMM right transpose"
          write(*,*) abs(sum(Mat-MatR))/M/M, " is larger then 1E-14!"
          stop
        endif
        
        !test right-mult op='c'
        DO I=1,M
        DO J=1,M
          Mat(I,J) = CMPLX(DBLE(i-j),DBLE(i+2*j),kind(0.d0))
        ENDDO
        ENDDO
        MatR(:,:)=Mat(:,:)
        DO I=1,M
        DO J=1,N
        Z = CMPLX(0.d0,0.d0,kind(0.d0))
        DO K=1,N
          Z = Z + Mat(I,P(K)) * conjg(A(j,k))
        ENDDO
        MatR(I,P(J)) = Z
        ENDDO
        ENDDO
        
        CALL ZSLGEMM('R','c', N, M, A, P, Mat)
        IF(abs(sum(Mat-MatR))/M/M>1E-14) then
          write(*,*) "ERROR in ZSLGEMM right conjg-transpose"
          write(*,*) abs(sum(Mat-MatR))/M/M, " is larger then 1E-14!"
          stop
        endif
        
!         enddo !O loop
        
        Deallocate(A,P)
        enddo ! N loop
        Deallocate(Mat,MatL,MatR)
        write(*,*) "SUCCESS"
        
END PROGRAM Main
