     Program Test

       Use MaxEnt_stoch_mod
       use Files_mod
       
       Implicit none
       
       Real (Kind=Kind(0.d0)), Dimension(:), allocatable :: XQMC, XTAU, Alpha_tot, xom, A, &
            &    ERROR, ARES, XQMC_ST
       Real (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: XCOV, Cov_st
       Real (Kind=Kind(0.d0)), External :: XKER, Back_trans_Aom, F_Fit
       Character (len=64) :: File1, File2
       Character (len=1)  :: Fermion_type
       Complex (Kind=Kind(0.d0)) :: Z
       Character (Len=64) :: command
       Logical            :: skip

       Integer ::  N_alpha, nc, nw, Ndis, nwp, nbins, NSweeps, NWarm, Ntau, Ngamma, nt, N_alpha_1, io_error
       Real (Kind=Kind(0.d0)) :: alpha_st, R, Xmom1, X,  X1, X2, om, omp, pi, Om_st, Om_en, err
       Real (Kind=Kind(0.d0)) :: Dom, delta, beta, dtau


       
       open(unit=30,file='paramSAC',status='old',action='read', iostat=io_error) 
       Open(unit=50,File='Info',Status="unknown")
       if (io_error.eq.0) then
          write(50,*) 'Reading in from paramSAC'
          read(30,*)  Ngamma, OM_st, OM_en, Ndis,  NBins, NSweeps, Nwarm
          read(30,*)  N_alpha, alpha_st, R
       else
          write(6,*) 'No file paramSAC! ' 
          stop
       endif
       close(30)
       

       open (unit=10,File="g_dat", status="unknown")
       Read(10,*) X, xmom1, X1
       if (X1 < 10.D-8) then
          Skip =.true.
          nc = 0
       else
          nc = 1
       endif
       do
          read(10,*,End=10)  X,X1,X2
          if ( abs(X1) - abs(X2)  < 0.d0 .or. X1 < 0.d0 ) exit
          nc = nc + 1
       enddo
10     continue
       Rewind(10)
       Write(6,*) Nc, Skip, X1, X2
       Ntau = NC
       Allocate ( XCOV(NTAU,NTAU), XQMC(NTAU), XQMC_ST(NTAU), XTAU(NTAU) )
       xcov = 0.d0
       if (Skip) read(10,*) X,X1,X2
       do nt = 1,ntau
          read(10,*) Xtau(nt), Xqmc(nt), err
          xcov(nt,nt) = err*err
       enddo
       close(10)
       xqmc_st = xqmc
       dtau = Xtau(2) - Xtau(1)
          
       Allocate (Alpha_tot(N_alpha) )
       do nt = 1,N_alpha
          alpha_tot(nt) = alpha_st*(R**(nt-1))
       enddo
       write(50,*) 'First Moment  ', Xmom1
       close(50)
          
          
       Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER, Back_Trans_Aom, Beta, &
            &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm) 
       
       
       N_alpha_1 = N_alpha - 10  
       File1 ="Aom_ps"
       file2 =File_i(File1,N_alpha_1)
       Write(50,*) file2
       Open(Unit=66,File=file2,status="unknown")
       Allocate (xom(Ndis), A(Ndis))
       do nw = 1,Ndis
          read(66,*) xom(nw), A(nw), x, x1, x2
       enddo
       close(66)
       Dom = xom(2) - xom(1)
       
       
       ! Compute the real frequency Green function.
       delta = 4.0*Dom
       Open (Unit=43,File="Green", Status="unknown", action="write")
       pi = acos(-1.d0)
       do nw = 1,Ndis
          Z = cmplx(0.d0,0.d0,Kind(0.d0))
          om = xom(nw)
          do nwp = 1,Ndis
             omp = xom(nwp)
             Z = Z + A(nwp)/cmplx( om -  omp, delta, kind(0.d0))
          enddo
          Z = Z * dom
          write(43,"('X'2x,F14.7,2x,F16.8,2x,F16.8)" ) xom(nw), dble(Z), -Aimag(Z)/pi
       enddo
       close(43)
       
       
       Open (Unit=70,File="data_out", status="unknown")
       do nt = 1,Ntau
          X = 0.d0
          do  nw = 1, Ndis
             X = X + Xker(Xtau(nt),xom(nw), beta)*A(nw)
          enddo
          X = X *dom
          write(70,2005) xtau(nt), xqmc_st(nt), sqrt(xcov(nt,nt)), X
       enddo
2005   format(F14.7,2x,F14.7,2x,F14.7,2x,F14.7)
       Close(70)
       
       Command = "rm dump*"
       Call System (Command)
       Command = "ls"
       Call System (Command)
       
          
     end Program Test
     
     
     
     Real (Kind=Kind(0.d0)) function XKER(tau,om, beta)
       
       Implicit None
       real (Kind=Kind(0.d0)) :: tau, om, pi, beta

       !pi = 3.1415927
       XKER = exp(-tau*om) !/ ( 1.d0 + exp(-Beta*om) ) ! /pi

     end function XKER

     Real (Kind=Kind(0.d0)) function Back_trans_Aom(Aom, om, beta)

       Implicit None
       real (Kind=Kind(0.d0)) ::  Aom, om, beta

       Back_trans_Aom = Aom 

     end function BACK_TRANS_AOM

     
     
