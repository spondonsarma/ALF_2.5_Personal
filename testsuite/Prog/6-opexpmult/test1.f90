! compile with
! gfortran -std=f2003  -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ main.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a

Program OPEXPMULTTEST

Use Operator_mod


        COMPLEX (KIND=KIND(0.D0)), DIMENSION(3,3) :: U
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(3,5) :: V
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(3) :: Z
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(5,5) :: matnew, matold
        Complex (KIND = KIND(0.D0)) :: tmp, lexp
        Integer, DIMENSION(3) :: P
        Integer :: i,j,n, opn, Ndim
        
        
        opn = 3
        Ndim = 5
        do i=1,Ndim
            do j=1,Ndim
                matnew(i,j) = CMPLX(i, j, kind(0.D0))
                matold(i,j) = CMPLX(i, j, kind(0.D0))
            enddo
        enddo
        
        do i = 1, opn
        do j = 1, Ndim
        V(i,j) = CMPLX(i, j, kind(0.D0))
        enddo
        P(i) = i
        Z(i) = exp(CMPLX(i, j, kind(0.D0)))
        enddo
        
        do i = 1, opn
        do j = 1, opn
        U(i,j) = CMPLX(i, j, kind(0.D0))
        enddo
        enddo
        
        call opexpmult(V, U, P, matnew, Z, opn, Ndim)
        
! check against old version
    do n = 1, opn
        lexp = Z(n)
        DO I = 1, Ndim
!             Mat(I,Op%P(n))  =  ExpMOp(n) * zdotu(Op%N,Op%U(1,n),1,VH(1,I),1) 
            tmp=cmplx(0.d0, 0.d0, kind(0.D0))
            Do m = 1, opn
                tmp = tmp + V(m,I) * U(m,n)
            Enddo
            Matold(I, P(n))  =  lexp * tmp
        enddo
    Enddo
    
    do i = 1,Ndim
        do j = 1,Ndim
        if (matold(i, j) .ne. matnew(i, j)) then
        STOP 2
        endif
        enddo
    enddo

end Program OPEXPMULTTEST
