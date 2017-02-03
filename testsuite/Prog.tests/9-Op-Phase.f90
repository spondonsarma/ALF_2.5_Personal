! compile with
!gfortran -std=f2003  -I ../../Libraries/Modules/ -I ../../Prog_8/  -L ../../Libraries/Modules/ 9-Op-Phase.f90 ../../Prog_8/Operator.o ../../Libraries/Modules/modules_90.a -llapack -lblas

Program TESTOP_PHASE

Use Operator_mod
implicit none

        Real (KIND = KIND(0.D0)) :: spin
        Complex(Kind = Kind(0.D0)) :: Phasenew, Phaseold, diff
        Integer :: i,n, N_SUN, NSigma(3,3), nf,nt
        Type(Operator) :: Op(3,3)

        call op_seths()
        do nf = 1, 3
            do nt = 1, 3
            N_SUN = 3
            Call Op_make(Op(nf,nt), 3)

            do i = 1, Op(nf, nt)%N
                Op(nf, nt)%E(i) = 2*i-3
            enddo
            Op(nf, nt)%N_non_zero = 2
            Op(nf, nt)%g = 2
            Op(nf,nt)%type = 1 + Mod(nf + nt, 2)
            Op(nf,nt)%alpha = cmplx(dble(nf), dble(nt), kind(0.D0))
            spin =-1.0
            nsigma(nf, nt) = (2*Mod(nf, 2) - 1)*(1 + Mod(nf+nt, 2))
            enddo
        enddo
        Phasenew = 1.D0
        Phaseold = 1.D0

        Call Op_Phase(Phasenew, Op, NSigma, N_SUN)

! check against old version
    do nf = 1, Size(Op,2)
       do n = 1, size(Op,1)
          do nt = 1, size(nsigma,2)
             Phaseold=Phaseold*exp(cmplx(0.d0,Aimag(Op(n,nf)%g*Op(n,nf)%alpha) * Phi(nsigma(n,nt),Op(n,nf)%type),kind(0.D0)))
          enddo
       enddo
    enddo
    Phaseold = Phaseold**N_SUN

    diff = Phaseold - Phasenew

    if (abs(DBLE(diff)) > MAX(ABS(DBLE(Phaseold)), ABS(DBLE(Phasenew)))*1D-14) then
    write (*,*) "ERROR in real", Phaseold, Phasenew, DBLE(diff), MAX(ABS(DBLE(Phaseold)), ABS(DBLE(Phasenew)))*1D-14
    STOP 2
    endif
    if (abs(AIMAG(diff)) > MAX(ABS(AIMAG(Phaseold)), ABS(AIMAG(Phasenew)))*2D-14) then
    write (*,*) "ERROR in imag", Phaseold, Phasenew, diff, MAX(ABS(AIMAG(Phaseold)), ABS(AIMAG(Phasenew)))*2D-14
    STOP 3
    endif
    write (*,*) "success"
end Program TESTOP_PHASE
