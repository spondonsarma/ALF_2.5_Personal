! compile with
! gfortran -Wall -std=f2003 -I ../../../Prog_8/  -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ main.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a


Program OPMULTTEST

Use Operator_mod
implicit none

        Complex (Kind=Kind(0.D0)) :: Z, Z1, Zre, Zim
        Real (KIND = KIND(0.D0)) :: spin
        Complex (Kind=Kind(0.D0)), Dimension(:, :), allocatable :: VH, matnew, matold
        Integer :: i, n, m, j, ndim, N_type, opn
        Type(Operator) :: Op
        
! setup some test data
        Ndim = 3
        
        do opn = 1, 3
        do N_type = 1, 2
        allocate(VH(opn, Ndim), matold(Ndim, Ndim), matnew(Ndim, Ndim))
        call op_seths()
        Call Op_make(Op, opn)
        
        do i = 1, Op%N
            Op%E(i) = 2*i-Ndim
            Op%P(i) = i
            do n = 1,Op%N
            Op%U(i,n) = CMPLX(n, i, kind(0.D0))
            enddo
        enddo
        Op%N_non_zero = max(1, opn - 1)
        Op%g = 2.D0
        spin =-1.0
        
        do i = 1,Ndim
        do n = 1,Ndim
        matnew(i,n) = CMPLX(i,n, kind(0.D0))
        matold(i,n) = CMPLX(i,n, kind(0.D0))
        enddo
        enddo
        
        Call Op_Wrapup(matnew, Op, spin, Ndim, N_type)

! check against old version from Operator_FFA.f90

   If (N_type == 1) then
       VH = 0.d0
       do n = 1,Op%N
          Z=cmplx(1.d0,0.d0, kind(0.D0))
          If ( n <= OP%N_non_Zero) Z = exp(-Op%g*cmplx(Op%E(n)*spin,0.d0, kind(0.D0))) 
          Do m = 1,Op%N
             Z1 = Op%U(m,n) * Z
             DO I = 1,Ndim
                VH(I,n)  = VH(I,n) + Matold(I,Op%P(m)) * Z1
             Enddo
          enddo
       Enddo
       Do n = 1,Op%N
          Do I = 1,Ndim
             Matold(I,Op%P(n)) =   VH(I,n) 
          Enddo
       Enddo

       VH = 0.d0
       do n = 1,Op%N
          Z=cmplx(1.d0,0.d0, kind(0.D0))
          If ( n <= OP%N_non_Zero) Z = exp(Op%g*cmplx(Op%E(n)*spin,0.d0, kind(0.D0)))
          Do m = 1,Op%N
             Z1 = Z * conjg(Op%U(m,n))
             DO I = 1,Ndim
                VH(I,n)  = VH(I,n) + Z1* Matold(Op%P(m),I) 
             Enddo
          enddo
       enddo
       Do n = 1,Op%N
          Do I = 1,Ndim
             Matold(Op%P(n),I) =   VH(I,n) 
          Enddo
       Enddo
    elseif (N_Type == 2) then
       VH = cmplx(0.d0,0.d0)
       do n = 1,Op%N
          Do m = 1,Op%N
             Z1 =  conjg(Op%U(n,m)) 
             DO I = 1,Ndim
                VH(I,n)  = VH(I,n) + Matold(I,Op%P(m)) * Z1
             Enddo
          enddo
       Enddo
       Do n = 1,Op%N
          Do I = 1,Ndim
             Matold(I,Op%P(n)) =   VH(I,n) 
          Enddo
       Enddo
       
       VH = cmplx(0.d0,0.d0)
       do n = 1,Op%N
          Do m = 1,Op%N
             Z1 =  Op%U(n,m)
             DO I = 1,Ndim
                VH(I,n)  = VH(I,n) + Z1* Matold(Op%P(m),I) 
             Enddo
          enddo
       enddo
       Do n = 1,Op%N
          Do I = 1,Ndim
             Matold(Op%P(n),I) =   VH(I,n) 
          Enddo
       Enddo
    endif

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

    deallocate(VH, matnew, matold)
enddo
enddo

end Program OPMULTTEST
