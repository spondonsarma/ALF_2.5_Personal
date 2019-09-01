     Program Test

       Use MaxEnt_stoch_mod

       Implicit  None
       !Implicit Real (Kind=Kind(0.d0)) (A-G,O-Z)
       !Implicit Integer (H-N)
       
       Real (Kind=Kind(0.d0)), Dimension(:)  , allocatable :: XQMC, XTAU, Alpha_tot, om_bf, alp_bf, xom, A
       Real (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: XCOV, Xn
       Real (Kind=Kind(0.d0))                              :: X_moments(2), Xerr_moments(2)
       Real (Kind=Kind(0.d0)), External                    :: XKER_ph, Back_trans_ph, XKER_pp, Back_trans_pp
       Character (Len=64)                                  :: command, File1, File2
       Complex (Kind=Kind(0.d0))                           :: Z
       

       Integer                :: Ngamma, Ndis,  NBins, NSweeps, Nwarm, N_alpha, L_cov
       Real (Kind=Kind(0.d0)) :: OM_st, OM_en,  alpha_st, R
       Logical                :: Checkpoint
       Character (Len=2)      :: Channel
       
       Integer                :: nt, nt1, io_error, n,nw, nwp, ntau, N_alpha_1, i
       Real (Kind=Kind(0.d0)) :: dtau, pi, xmom1, x,x1,x2, tau, omp, om, Beta,err, delta, Dom

       NAMELIST /VAR_Max_Stoch/ Ngamma, Ndis,  NBins, NSweeps, Nwarm, N_alpha, L_cov, &
            &                   OM_st, OM_en,  alpha_st, R,  Checkpoint, Channel

       open(unit=30,file='paramSAC',status='old',action='read', iostat=io_error) 
       if (io_error.eq.0) then
          READ(30,NML=VAR_Max_Stoch)
       else
          write(6,*) 'No file paramSAC ' 
          stop
       endif

       Open(unit=50,File='Info_MaxEnt',Status="unknown")
       write(50,*) 'Channel      :: ', Channel
       If (Channel == "ph" )  then
          Write(50,*)  'Om_start is set to zero. PH channel corresponds to symmetric data '
          Om_st = 0.d0
       endif
       Write(50,*) 'Covariance         :: ', L_cov
       Write(50,*) 'Checkpoint         :: ', Checkpoint
       Write(50,*) 'Om_st, Om_en       :: ', Om_st, Om_en
       Write(50,*) 'Delta Om           :: ', (Om_en - Om_st)/real(Ndis,kind(0.d0))
       Write(50,*) 'Bins, Sweeps, Warm :: ', NBins, NSweeps, Nwarm
       If (N_alpha <= 10 ) then
          Write(50,*) 'Not enough temperatures: N_alpha has to be bigger than 10'
          Stop
       Endif
       Write(50,*) 'N_Alpha, Alpha_st,R:: ', N_alpha, alpha_st, R
       

       open (unit=10,File="g_dat", status="unknown") 
       read(10,*)  ntau, n, Beta
       Allocate ( XCOV(NTAU,NTAU), XQMC(NTAU),XTAU(NTAU) )
       XCOV  = 0.d0
       Do nt = 1,NTAU
          read(10,*)  xtau(nt), xqmc(nt), err 
          xcov(nt,nt) = err*err
       Enddo
       if (L_cov.eq.1) then
          do nt = 1,ntau
             do nt1 = 1,ntau
                read(10,*) xcov(nt,nt1)
             enddo
          enddo
       endif
       close(10)
       dtau = Xtau(2) - Xtau(1)

       pi = acos(-1.d0)
       Select Case (Channel)
       Case ("PH")
          xmom1 = pi * xqmc(1) 
       Case ("PP")
          xmom1 = 2.d0* pi * xqmc(1)
       Case ("P")
          xmom1 =  pi 
       Case default 
          Write(6,*) "Channel not yet implemented"
          Stop
       end Select

       Allocate (Alpha_tot(N_alpha) )
       do nt = 1,N_alpha
          alpha_tot(nt) = alpha_st*(R**(nt-1))
       enddo
       write(50,*) 'First Moment, Beta  ', Xmom1, Beta
       close(50)

       Select Case (Channel)
       Case ("PH")
          If (L_cov == 1 ) then
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_ph, Back_Trans_ph, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm, L_Cov)
          else
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_ph, Back_Trans_ph, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm)
          endif
          ! Beware: Xqmc and cov are modified in the MaxEnt_stoch call.
       Case ("PP")
          If (L_cov == 1 ) then
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_pp, Back_Trans_pp, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm, L_Cov)
          else
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_pp, Back_Trans_pp, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm)
          endif
          ! Beware: Xqmc and cov are modified in the MaxEnt_stoch call.
       Case default 
          Write(6,*) "Channel not yet implemented"
          Stop
       end Select
       
       If ( .not.  Checkpoint ) then
          Command = "rm dump*"
          Call System (Command)
          Command = "ls"
          Call System (Command)
       endif
       
       open (unit=10,File="g_dat", status="unknown") 
       read(10,*)  ntau, n, beta
       XCOV  = 0.d0
       Do nt = 1,NTAU
          read(10,*)  xtau(nt), xqmc(nt), err 
          xcov(nt,nt) = err*err
       Enddo
       close(10)
       pi = acos(-1.d0)
       xmom1 = pi * xqmc(1) 

       Open (Unit = 10,File="Best_fit", Status ="unknown") 
       Allocate (om_bf(Ngamma), alp_bf(Ngamma) )
       DO i = 1, Ngamma
          read(10,*)  om_bf(i), alp_bf(i)
       Enddo
       close(10) 

       Open (Unit = 11,File="Data_out", Status ="unknown") 
       do nt = 1,Ntau
          X = 0.d0
          tau = xtau(nt)
          do i = 1,Ngamma
             X = X + alp_bf(i)*Xker_ph(tau,om_bf(i), beta)
          enddo
          Write(11,"(F14.7,2x,F14.7,2x,F14.7,2x,F14.7)")  xtau(nt), xqmc(nt),  sqrt(xcov(nt,nt)), xmom1*X
       enddo
       close(11)

       N_alpha_1 = N_alpha - 10  
       File1 ="Aom_ps"
       file2 =File_i(File1,N_alpha_1)
       !Write(50,*) file2
       Open(Unit=66,File=file2,status="unknown")
       Allocate (xom(Ndis), A(Ndis))
       do nw = 1,Ndis
          read(66,*) xom(nw), A(nw), x, x1, x2
       enddo
       close(66)
       Dom = xom(2) - xom(1)
       
       ! Compute the real frequency Green function.
       delta = Dom
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
       
       
     end Program Test
     


     Real (Kind=Kind(0.d0)) function XKER_ph(tau,om, beta)

       Implicit None
       real (Kind=Kind(0.d0)) :: tau, om, pi, beta

       pi = 3.1415927

       XKER_ph = (exp(-tau*om) + exp(-( beta - tau )*om ) )/( pi*(1.d0 + exp( - beta * om ) ) )

     end function XKER_ph

     Real (Kind=Kind(0.d0)) function XKER_pp(tau,om, beta)

       Implicit None
       real (Kind=Kind(0.d0)) :: tau, om, pi, beta

       pi = 3.1415927

       XKER_pp = exp(-tau*om) / ( pi*(1.d0 + exp( - beta * om ) ) )

     end function XKER_pp

     Real (Kind=Kind(0.d0)) function XKER_p(tau,om, beta)

       Implicit None
       real (Kind=Kind(0.d0)) :: tau, om, pi, beta

       pi = 3.1415927

       XKER_p  = exp(-tau*om) / ( pi*(1.d0 + exp( - beta * om ) ) )

     end function XKER_p

     
     Real (Kind=Kind(0.d0)) function Back_trans_ph(Aom, om, beta)

       Implicit None
       real (Kind=Kind(0.d0)) ::  Aom, om, beta

       Back_trans_ph = Aom/(1.d0 + exp(-beta*om) )   
       ! This gives S(q,om) = chi(q,om)/(1 - e^(-beta om))

     end function BACK_TRANS_PH

     Real (Kind=Kind(0.d0)) function Back_trans_pp(Aom, om, beta)
            
       Implicit None
       real (Kind=Kind(0.d0)) ::  Aom, om, beta

       real (Kind=Kind(0.d0)) :: Zero

       Zero = 1.D-8

       if ( abs(om) < zero ) then
          Back_trans_pp = beta * Aom/2.d0
       else
          Back_trans_pp = Aom * (1.d0 + exp(-beta*om) ) / (om *( 1.d0 + exp(-beta*om) ) )
       endif
       ! This gives  = chi(q,om)/omega

     end function BACK_TRANS_PP

     Real (Kind=Kind(0.d0)) function Back_trans_p(Aom, om, beta)
            
       Implicit None
       real (Kind=Kind(0.d0)) ::  Aom, om, beta

       Back_trans_p =  Aom
       
     end function BACK_TRANS_P

     
