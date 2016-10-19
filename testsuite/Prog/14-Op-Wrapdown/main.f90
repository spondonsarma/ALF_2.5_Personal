! compile with
! gfortran -Wall -std=f2003 -I ../../../Prog_8/  -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ main.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a


Program OPMULTTEST

     Use Operator_mod
      Implicit None
!
      Interface
         Subroutine Op_WrapdoFFA (Mat, Op, spin, Ndim, N_Type)
            Use Operator_mod
            Type (Operator), Intent (In) :: Op
            Complex (Kind=8), Intent (Inout) :: Mat (Ndim, Ndim)
            Real (Kind=8), Intent (In) :: spin
            Integer, Intent (In) :: N_Type, Ndim
         End Subroutine
      End Interface

        Complex (Kind=Kind(0.D0)) :: Z, Z1, Zre, Zim, tmp, exphere
        Real (KIND = KIND(0.D0)) :: spin
        Complex (Kind=Kind(0.D0)), Dimension(:, :), allocatable :: VH, matnew, matold
        Complex (Kind=Kind(0.D0)), Dimension(:), allocatable :: Expop, ExpMop
        Integer :: i, n, m, j, ndim, N_type, opn
        Type(Operator) :: Op
        
! setup some test data
        Ndim = 30
        
      Do opn = 1, 4
         Do N_Type = 1, 2
            Allocate (VH(opn, Ndim), matold(Ndim, Ndim), matnew(Ndim, &
           & Ndim), Expop(opn), ExpMop(opn))
            Call Op_seths ()
            Call Op_make (Op, opn)
!
            Do i = 1, Op%n
               Op%P (i) = i
               Do n = 1, Op%n
                  Op%O (i, n) = CMPLX (n+i, 0.D0, kind(0.D0))
               End Do
            End Do
!
            Op%g = 2.D0
            Op%alpha = 0.D0
            Call Op_set (Op)
!
            spin = - 1.0
!
            Do i = 1, Ndim
               Do n = 1, Ndim
                  matnew (i, n) = CMPLX (i, n, kind(0.D0))
                  matold (i, n) = CMPLX (i, n, kind(0.D0))
               End Do
            End Do
        
        Call Op_Wrapdo(matnew, Op, spin, Ndim, N_type)

! check against old version from Operator_FFA.f90
call Op_WrapdoFFA(matold, Op, spin, Ndim, N_type)

    do i=1,3
    do j=1,3
    Zre = real(matnew(i,j)-matold(i,j))
    Zim = aimag(matnew(i,j)-matold(i,j))
    if (Abs(Zre) > MAX(ABS(real(matnew(i,j))), ABS(real(matold(i,j))) )*1D-15) then
    write (*,*) "opn: ", opn, "N_type", N_type
    write (*,*) "error in real part", real(matnew(i,j)), real(matold(i,j))
    STOP 2
    endif
    if (Abs(Zim) > MAX(ABS(aimag(matnew(i,j))), ABS(aimag(matold(i,j))) )*1D-15) then
    write (*,*) "error in imag part", aimag(matnew(i,j)), aimag(matold(i,j))
    STOP 3
    endif
    enddo
    enddo

    deallocate(VH, matnew, matold, expop, ExpMop)
enddo
enddo

end Program OPMULTTEST

  Subroutine Op_WrapdoFFA(Mat,Op,spin,Ndim,N_Type)

  Use Operator_mod
    Implicit none 

    Integer, intent(in) :: Ndim
    Type (Operator) , INTENT(IN )   :: Op
    Complex (Kind=8), INTENT(INOUT) :: Mat (Ndim,Ndim)
    Real    (Kind=8), INTENT(IN )   :: spin
    Integer, INTENT(IN) :: N_Type

    ! Local 
    Complex (Kind=8) :: VH(Ndim,Op%N), Z, Z1
    Integer :: n, i, m, m1
    
    !!!!! N_Type == 1
    !    Op%U*exp(-Op%g*spin*Op%E)*Mat*exp(Op%g*spin*Op%E)*(Op%U^{dagger})
    !    
    !!!!!
    !!!!! N_Type == 2
    !    (Op%U^{dagger}) * Mat * Op%U
    !!!!!
    If (N_type == 1) then
       VH = cmplx(0.d0,0.d0)
       Do m = 1,Op%N
          Z = cmplx(1.d0,0.d0, kind(0.D0))
          If ( m <= OP%N_non_Zero) Z = exp(Op%g*cmplx(Op%E(m)*spin,0.d0, kind(0.D0))) 
          do n = 1,Op%N
             Z1 = Z * conjg(Op%U(n,m))
             DO I = 1,Ndim
                VH(I,n)  = VH(I,n) + Mat(I,Op%P(m)) * Z1
             Enddo
          enddo
       Enddo
       Do n = 1,Op%N
          Do I = 1,Ndim
             Mat(I,Op%P(n)) =   VH(I,n) 
          Enddo
       Enddo

       VH = cmplx(0.d0,0.d0)
       Do m = 1,Op%N
          Z = cmplx(1.d0,0.d0, kind(0.D0))
          If ( m <= OP%N_non_Zero) Z = exp(-Op%g*cmplx(Op%E(m)*spin,0.d0, kind(0.D0)))
          do n = 1,Op%N
             Z1 = Z * Op%U(n,m)
             DO I = 1,Ndim
                VH(I,n)  = VH(I,n) + Z1* Mat(Op%P(m),I) 
             Enddo
          enddo
       enddo
       Do n = 1,Op%N
          Do I = 1,Ndim
             Mat(Op%P(n),I) =   VH(I,n) 
          Enddo
       Enddo
    elseif (N_Type == 2) then
       VH = cmplx(0.d0,0.d0, kind(0.D0))
       do n = 1,Op%N
          Do m = 1,Op%N
             Z1 =  Op%U(m,n)
             DO I = 1,Ndim
                VH(I,n)  = VH(I,n) + Mat(I,Op%P(m)) * Z1
             Enddo
          enddo
       Enddo
       Do n = 1,Op%N
          Do I = 1,Ndim
             Mat(I,Op%P(n)) =   VH(I,n) 
          Enddo
       Enddo
       
       VH = cmplx(0.d0,0.d0, kind(0.D0))
       do n = 1,Op%N
          Do m = 1,Op%N
             Z1 =  conjg(Op%U(m,n))
             DO I = 1,Ndim
                VH(I,n)  = VH(I,n) + Z1* Mat(Op%P(m),I) 
             Enddo
          enddo
       enddo
       Do n = 1,Op%N
          Do I = 1,Ndim
             Mat(Op%P(n),I) =   VH(I,n) 
          Enddo
       Enddo
    endif
    
  end Subroutine Op_WrapdoFFA

