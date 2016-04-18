      SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
C
      INTEGER N,NM,IERR,MATZ
      DOUBLE PRECISION A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     OF A REAL SYMMETRIC MATRIX.
C
C     ON INPUT
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT.
C
C        N  IS THE ORDER OF THE MATRIX  A.
C
C        A  CONTAINS THE REAL SYMMETRIC MATRIX.
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT
C
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.
C
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT
C           AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO.
C
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  TRED1(NM,N,A,W,FV1,FV2)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
      END
