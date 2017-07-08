! compile with
! gfortran -Wall -std=f2003 -I ../../../Prog/  -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ 21-Op-make.f90 ../../../Prog/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas
!
Program TESTOPMAKE
!
      Use Operator_mod
      Implicit None
!
      Interface
        subroutine Op_exp_FFA(g,Op,Mat)
            Use Operator_mod
            Type (Operator), Intent(IN)  :: Op
            Complex (Kind=8), Dimension(:,:), INTENT(OUT) :: Mat
            Complex (Kind=8), INTENT(IN) :: g
         End Subroutine
      End Interface
!
      Complex (Kind=Kind(0.D0)) :: Zre, Zim
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: VH, &
     & matnew, matold
      Complex (Kind=Kind(0.D0)), Dimension (:), Allocatable :: Expop, &
     & ExpMop
      Integer :: i, n, j, Ndim, opn
      Type (Operator) :: Op
!
! setup some test data
      Ndim = 5
!
      Do opn = 1, 8
      Do Ndim = opn+1, 15
            Allocate (VH(opn, Ndim), matold(Ndim, Ndim), matnew(Ndim, &
           & Ndim), Expop(opn), ExpMop(opn))
            Call Op_seths ()
            Call Op_make (Op, opn)
!
            Do i = 1, Op%n
               Op%P (i) = i
               Do n = 1, Op%n
                  Op%O (i, n) = CMPLX (n+i, i-n, kind(0.D0))
               End Do
            End Do
!
            Op%g = 2.D0
            Op%alpha = 0.D0
            Call Op_set (Op)
! check against old version from Operator_FFA.f90
!
            Call Op_exp_FFA (Op%g, Op, matold)
!
            Call Op_exp (Op%g, Op, matnew)
!
!
            Do i = 1, Ndim
               Do j = 1, Ndim
                  Zre = real (matnew(i, j)-matold(i, j))
                  Zim = aimag (matnew(i, j)-matold(i, j))
                  If (Abs(Zre) > Max(Abs(real(matnew(i, j))), &
                 & Abs(real(matold(i, j))))*1D-14) Then
                     Write (*,*) "opn: ", opn
                     Write (*,*) "ERROR in real part", real (matnew(i, &
                    & j)), real (matold(i, j))
                     Stop 2
                  End If
                  If (Abs(Zim) > Max(Abs(aimag(matnew(i, j))), &
                 & Abs(aimag(matold(i, j))))*1D-14) Then
                     Write (*,*) "opn: ", opn
                     Write (*,*) "ERROR in imag part", aimag (matnew(i, &
                    & j)), aimag (matold(i, j))
                     Stop 3
                  End If
               End Do
            End Do
!
            Deallocate (VH, matnew, matold, Expop, ExpMop)
            call Op_clear(Op, opn)
         End Do
      End Do

End Program TESTOPMAKE
!
  subroutine Op_exp_FFA(g,Op,Mat)
    Use Operator_mod
    Implicit none 
    Type (Operator), Intent(IN)  :: Op
    Complex (Kind=8), Dimension(:,:), INTENT(OUT) :: Mat
    Complex (Kind=8), INTENT(IN) :: g
    Complex (Kind=8) :: Z, Z1

    Integer :: n, i,j
      
    Mat = cmplx(0.d0,0.d0)
    Do n = 1,Op%N
       Z = exp(g*cmplx(Op%E(n),0.d0))
       do J = 1,Op%N
          Z1 = Z*conjg(Op%U(J,n))
          Do I = 1,Op%N
             Mat(I,J) = Mat(I,J) + Op%U(I,n)*Z1
          enddo
       enddo
    enddo
  end subroutine Op_exp_FFA
