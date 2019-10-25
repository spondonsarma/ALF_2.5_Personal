!  Copyright (C) 2019 The ALF project
! 
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
!     You should have received a copy of the GNU General Public License
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
!     
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!     
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de .
!       
!     - We require the preservation of the above copyright notice and this license in all original files.
!     
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
! 
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.

   module ana_mod 
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Collection of routines for postprocessing the Monte Carlo bins
!
!--------------------------------------------------------------------
      Use Errors
      Use MyMats
      Use Matrix
      Use Lattices_v3

   contains

   Subroutine read_vec(file, sgn, bins)
!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief
!> Reads in bins of scalar observables from file
!> 
!> @param [IN] file Character(len=64)
!> \verbatim
!>  Name of file that gets read in
!> \endverbatim
!> @param [OUT] sgn Real(:)
!> \verbatimam
!>  Sign of bins 
!> \endverbatim
!> @param [OUT] bins Complex(:,:)
!> \verbatim
!>  Monte Carlo bins
!> \endverbatim
!-------------------------------------------------------------------
      Implicit none
      Character (len=64), intent(in) :: file
      Real    (Kind=Kind(0.d0)), allocatable, intent(out) :: sgn(:)
      Complex (Kind=Kind(0.d0)), pointer, intent(out) :: bins(:,:)
      
      Integer :: N, N1, I, Nobs, Nbins
      Real    (Kind=Kind(0.d0)) :: X
      Complex (Kind=Kind(0.d0)), Allocatable  :: Tmp(:)
      
      
      Open (Unit=10, File=file, status="unknown")
      Read(10,*) NOBS
      NOBS = NOBS-1
      allocate (Tmp(NOBS) )
      rewind(10)
      Nbins = 0
      do
         read(10,*,End=10) N,  (Tmp(I), I=1,size(Tmp,1)), X
         Nbins = Nbins + 1
      enddo
10    continue
      Close(10) 
      
      ALLOCATE(bins(NOBS,Nbins))
      ALLOCATE(sgn(Nbins))

      OPEN (UNIT=20, FILE=file, STATUS='old')
      DO N = 1,Nbins
         READ(20,*) N1, (Tmp(I), I=1,size(Tmp,1)),X 
         
         bins(:,N) = Tmp(:)
         sgn(N)  = X
      ENDDO
      CLOSE(20)
      deallocate(tmp)
   End Subroutine read_vec

!==============================================================================

   Subroutine read_latt(file, sgn, bins, bins0, Latt, dtau)
!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief
!> Reads in bins of lattice-type observables (both equal time and timedisplaced) from file
!> 
!> @param [IN] file Character(len=64)
!> \verbatim
!>  Name of file that gets read in
!> \endverbatim
!> @param [OUT] sgn Real(:)
!> \verbatimam
!>  Sign of bins 
!> \endverbatim
!> @param [OUT] bins Complex(:,:,:,:,:)
!> \verbatim
!>  Monte Carlo bins of correlation
!> @param [OUT] bins0 Complex(:,:)
!> \verbatim
!>  Monte Carlo bins of background
!> @param [OUT] Latt Type(Lattice)
!> \verbatim
!>  Bravais lattice
!> \endverbatim
!> @param [OUT] dtau Real
!> \verbatim
!>  Size of imaginary time step
!> \endverbatim
!-------------------------------------------------------------------
      Implicit none
      Character (len=64), intent(in) :: file
      Real    (Kind=Kind(0.d0)), allocatable, intent(out) :: sgn(:)
      Complex (Kind=Kind(0.d0)), pointer, intent(out) :: bins(:,:,:,:,:)
      Complex (Kind=Kind(0.d0)), pointer, intent(out) :: bins0(:,:)
      Type (Lattice), intent(out)    :: Latt
      Real    (Kind=Kind(0.d0)), intent(out) :: dtau

      Integer :: no, no1, n, nt, nb, Ntau, Nunit, Nbins, Norb!, LT
      Real    (Kind=Kind(0.d0)):: X, Y !, dt
      Real    (Kind=Kind(0.d0)), allocatable :: Xk_p(:,:)
      Real    (Kind=Kind(0.d0))              :: x_p(2)
      Complex (Kind=Kind(0.d0)) :: Z
      Real    (Kind=Kind(0.d0)), allocatable :: a1_p(:), a2_p(:), L1_p(:), L2_p(:)
      
      Integer :: ndim 
      
      Integer             :: L1, L2, N_SUN, Model_vers, ierr
      Character (len=64)  :: Model, Lattice_type
      Logical             :: Checkerboard, Symm
      NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model, N_SUN, Checkerboard, Symm, Model_vers
      
      
      ! Read in lattice
      open ( Unit=10, File=file, status="unknown" )
      read(10,*,ERR=100) X, Norb, Nunit, Ntau, dtau, ndim
      rewind(10)
      allocate( a1_p(ndim), a2_p(ndim), L1_p(ndim), L2_p(ndim) )
      read(10,*) X, Norb, Nunit, Ntau, dtau, ndim, L1_p, L2_p, a1_p, a2_p
      close(10)
      goto 120 
      
  100 continue
      Ntau = 1
      dtau = -1.d0
      rewind(10)
      Read(10,*,ERR=110) X, Norb, Nunit, Ntau, dtau
  110 continue
      rewind(10)
      Read(10,*) X, Norb, Nunit
      Close(10)
      
      allocate( a1_p(2), a2_p(2), L1_p(2), L2_p(2) )
      OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read')
      READ(5,NML=VAR_lattice)
      CLOSE(5)
      If ( Lattice_type =="BipartiteSquare" ) then
         a1_p(1) =  1.D0/sqrt(2.D0)  ; a1_p(2) =  1.D0/sqrt(2.D0)
         a2_p(1) =  1.D0/sqrt(2.D0)  ; a2_p(2) = -1.D0/sqrt(2.D0)
         L1_p    =  dble(L1)*a1_p
         L2_p    =  dble(L2)*a2_p
      elseif ( Lattice_type =="Square" ) then
         a1_p(1) =  1.0  ; a1_p(2) =  0.d0
         a2_p(1) =  0.0  ; a2_p(2) =  1.d0
         L1_p    =  dble(L1)*a1_p
         L2_p    =  dble(L2)*a2_p
      elseif ( Lattice_type=="Honeycomb" ) then
         a1_p(1) =  1.d0   ; a1_p(2) =  0.d0
         a2_p(1) =  0.5d0  ; a2_p(2) =  sqrt(3.d0)/2.d0
         L1_p    =  dble(L1) * a1_p
         L2_p    =  dble(L2) * a2_p
      elseif ( Lattice_type == "Pi_Flux" ) then 
         a1_p(1) =  1.D0   ; a1_p(2) =   1.d0
         a2_p(1) =  1.D0   ; a2_p(2) =  -1.d0
         !del_p   =  (a2_p - 0.5*a1_p ) * 2.0/3.0
         L1_p    =  dble(L1) * (a1_p - a2_p)/2.d0
         L2_p    =  dble(L2) * (a1_p + a2_p)/2.d0
      else
         Stop "Lattice not yet implemented!"
      endif
      
  120 continue
      Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
      deallocate( a1_p, a2_p, L1_p, L2_p )

      ! Determine the number of bins. 
      Open ( Unit=10, File=file, status="unknown" )
      nbins = 0
      do
         Read(10,*,End=10)
         Do no = 1,Norb
            Read(10,*)
         enddo
         do n = 1,Nunit
            Read(10,*)
            do nt = 1,Ntau
               do no = 1,Norb
                  do no1 = 1,Norb
                     read(10,*)
                  enddo
               enddo
            enddo
         enddo
         nbins = nbins + 1
      enddo
   10 continue
      Close(10)
      
      ! Allocate  space
      Allocate ( bins(Nunit,Ntau,Norb,Norb,Nbins), sgn(Nbins), Xk_p(2,Nunit), bins0(Norb,Nbins))
      
      Open ( Unit=10, File=file, status="unknown" ) 
      do nb = 1, nbins
         Read(10,*,End=10) sgn(nb)
         Do no = 1,Norb
            Read(10,*) bins0(no,nb)
         Enddo
         do n = 1,Nunit
            Read(10,*) Xk_p(1,n), Xk_p(2,n)
            do nt = 1,Ntau
               do no = 1,norb
                  do no1 = 1,Norb
                     read(10,*) bins(n,nt,no,no1,nb)
                  enddo
               enddo
            enddo
         enddo
      enddo
      close(10)
      
      do n = 1,Nunit
         x_p = dble(Latt%listk(n,1))*Latt%b1_p + dble(Latt%listk(n,2))*Latt%b2_p  
         x   = (x_p(1)-Xk_p(1,n))**2 + (x_p(2)-Xk_p(2,n))**2
         if ( x > 0.00001 ) then
            print*, "Error in read_latt: momenta do not fit", x, x_p, Xk_p(1,n)
            stop
         endif
      enddo
   
   End Subroutine read_latt

!==============================================================================

   Subroutine jack(func, data, N_skip, N_rebin, best, err)
!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief
!> Performs jackknife error Analysis
!> 
!-------------------------------------------------------------------
      Implicit none
      
      Complex (Kind=Kind(0.d0)), External                :: func
      Complex (Kind=Kind(0.d0)), intent(in)              :: data(:,:)
      Integer, intent(in)                                :: N_skip, N_rebin
      
      Complex (Kind=Kind(0.d0)), intent(out) :: best, err
      
      Complex (Kind=Kind(0.d0)), allocatable :: data_r(:,:), j(:), X(:), X2(:)
      Real (Kind=Kind(0.d0)), allocatable :: Rhelp(:)
      real (Kind=Kind(0.d0)) :: XM, XERR
      Integer :: Nobs, Nbins, Nbins_r
      Integer :: N_r, N1, N
  
  
      Nobs  = size(data,1)
      Nbins = size(data,2)
  
      ! Skipping and Rebinning 
      Nbins_r = (Nbins-N_skip)/N_rebin
      Allocate( data_r(Nobs,Nbins_r), X(Nobs) )
      N = N_skip
      Do N_r = 1, Nbins_r
         X(:) = cmplx(0.D0,0.d0,kind(0.d0))
         Do N1 = 1, N_rebin
            N = N + 1
            X(:) = X(:) + data(:,N)
         enddo
         data_r(:,N_r) = X(:)/dble(N_rebin)
      enddo
  
      ! Creating jackknife bins
      Allocate( j(Nbins_r), X2(Nobs) )
      X(:) = cmplx(0.D0,0.d0,kind(0.d0))
      do N = 1, Nbins_r
         X(:) = X(:) + data_r(:,N)
      enddo
      do N = 1, Nbins_r
         X2(:) = (X(:) - data_r(:,N))/dble(Nbins_r-1)
         j(N) = func(X2)
      enddo
  
      ! Calculate standard deviation of jackknife bins
      Allocate (Rhelp(Nbins_r))
      do N = 1, Nbins_r
         Rhelp(N) = dble(j(N))
      enddo
      call errcalc(Rhelp, xm, xerr)
      best =  cmplx(xm  , 0.d0, kind(0.D0))
      err  =  cmplx(xerr, 0.d0, kind(0.D0))

      do N = 1, Nbins_r
         Rhelp(N) = aimag(j(N))
      enddo
      call errcalc(Rhelp, xm, xerr)
      best =  best + cmplx( 0.d0, xm, kind(0.D0)   )
      err  =  err  + cmplx( 0.d0, xerr, kind(0.D0) )
      
      err = err * dble(Nbins_r)
      
      deallocate( data_r, X, j, X2, Rhelp )
  
   End Subroutine jack

!==============================================================================

   subroutine auto(func, data, N_skip, res)
!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief
!> Calculates autocorrelation
!-------------------------------------------------------------------
      implicit none
      
      complex (Kind=Kind(0.d0)), External   :: func
      complex (Kind=Kind(0.d0)), intent(in) :: data(:,:)
      Integer, intent(in)                   :: N_skip
      
      Real (Kind=Kind(0.d0)), intent(inout) :: res(:)
      
      Integer                :: N_obs, N_bins, N_auto, ntau, nt
      complex (Kind=Kind(0.d0)), allocatable :: data1(:), data2(:)
      complex (Kind=Kind(0.d0)) :: mean, X1, X2
      
      N_obs  = size(data,1)
      N_bins = size(data,2) - N_skip
      N_auto = size(res)
      allocate( data1(N_obs), data2(N_obs) )
      
      do ntau = 1, N_auto
         X1 = 0.0
         X2 = 0.0
         mean = 0.0
         do nt = 1, N_bins - ntau
            data1(:) = data(:,nt+N_skip)
            mean = mean + func(data1)
         enddo
         mean = mean / dble(N_bins - ntau)
         mean = func(mean)
         
         do nt = 1, N_bins - ntau
            data1(:) = data(:,nt+N_skip)
            data2(:) = data(:,nt+N_skip+ntau)
            X1 = X1 + (func(data1)-mean)*(func(data2)-mean) 
            X2 = X2 + (func(data1)-mean)*(func(data1)-mean)
         enddo
      !     X1 = X1 / dble(N_bins - ntau)
      !     X2 = X2 / dble(N_bins - ntau)
                  
         Res(ntau) = dble( X1/X2 )
      enddo
      
   
   end subroutine auto

!==============================================================================


   subroutine Cov_tau(file, PartHole)
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Analysis of imaginary time displaced  correlation functions.
!> 
!> @param [IN] file Character(len=64)
!> \verbatim
!>  Name of file that gets read in
!> \endverbatim
!> @param [IN] PartHole logical
!> \verbatimam
!>  If true, it is assumed that observable is symmetric in imaginary time around tau = beta/2
!> \endverbatim
!
!--------------------------------------------------------------------        

      Implicit none
      Character (len=64), intent(in) :: file
      Logical           , intent(in) :: PartHole
      
      Character (len=64) :: name_obs
      Real    (Kind=Kind(0.d0)), allocatable :: sgn(:)
      Complex (Kind=Kind(0.d0)), pointer :: Bins_raw(:,:,:,:,:), Bins0_raw(:,:)
      Type (Lattice)   :: Latt
      real    (Kind=Kind(0.d0)):: dtau
      Integer :: i

      i = len(trim(file)) - 4
      name_obs = file(:i)
      
      call read_latt(file, sgn, bins_raw, bins0_raw, Latt, dtau)
      call ana_tau(name_obs, sgn, bins_raw, bins0_raw, Latt, dtau, PartHole)

   end subroutine Cov_tau

!==============================================================================

   Subroutine ana_tau(name_obs, sgn, bins_raw, bins0_raw, Latt, dtau, PartHole)
      Implicit none
      Character (len=64), intent(in) :: name_obs
      Real    (Kind=Kind(0.d0)), allocatable, intent(in) :: sgn(:)
      Complex (Kind=Kind(0.d0)), pointer, intent(in) :: Bins_raw(:,:,:,:,:)
      Complex (Kind=Kind(0.d0)), pointer, intent(in) :: Bins0_raw(:,:)
      Type (Lattice), intent(in)    :: Latt
      Real    (Kind=Kind(0.d0)), intent(in) :: dtau
      Logical, intent(in) :: PartHole
      
      Character (len=64) :: File_out
      Real    (Kind=Kind(0.d0)), parameter :: Zero=1.D-8
      Integer :: N_skip, N_rebin, N_Cov, N_Back, N_auto
      Integer :: Nbins, Norb, LT, Lt_eff, Nunit
      Integer :: nb, no, no1, n, nt, nt1, ierr
      Complex (Kind=Kind(0.d0)) :: Z, Zmean, Zerr
      Real    (Kind=Kind(0.d0)), allocatable :: Phase(:)
      Complex (Kind=Kind(0.d0)), allocatable :: PhaseI(:)
      Complex (Kind=Kind(0.d0)), allocatable :: Bins(:,:,:)
      Real    (Kind=Kind(0.d0)), allocatable :: Xk_p(:,:)
      Complex (Kind=Kind(0.d0)), allocatable :: V_help(:,:)
      Complex (Kind=Kind(0.d0)), allocatable :: Bins_chi(:,:)
      Complex (Kind=Kind(0.d0)), allocatable :: Xmean(:), Xcov(:,:)
      
      NAMELIST /VAR_errors/   n_skip, N_rebin, N_Cov, N_Back, N_auto
      
      N_Back = 1
      N_auto = 0
      OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'unable to open <parameters>',ierr
         STOP
      END IF
      READ(5,NML=VAR_errors)
      CLOSE(5)
      
      Nbins = size(bins_raw,5)
      Norb  = size(bins_raw,3)
      LT    = size(bins_raw,2)
      Nunit = Latt%N
      
      Write(6,*) "# of bins: ", Nbins
      nbins  = Nbins - n_skip
      Write(6,*) "Effective # of bins: ", Nbins/N_rebin
      if(Nbins/N_rebin < 2) then
         stop "Effective # of bins smaller than 2. Analysis impossible!"
      endif

      if (PartHole .and. mod(Lt-1,2) == 0 ) then
         Lt_eff = (Lt -1 ) /2   + 1
      elseif (PartHole) then
         Lt_eff = Lt/2 
      else
         Lt_eff = Lt
      endif
      
      ! Allocate  space
      Allocate ( bins(Nunit,Lt_eff,Nbins), Phase(Nbins),PhaseI(Nbins), Xk_p(2,Nunit), &
            &     V_help(Lt_eff,Nbins))
      Allocate ( Bins_chi(Nunit,Nbins) )
      Allocate (Xmean(Lt_eff), Xcov(Lt_eff,Lt_eff))
      bins  = cmplx(0.d0,0.d0,Kind(0.d0))
   
      do n = 1,Nunit
         xk_p(:,n) = dble(Latt%listk(n,1))*Latt%b1_p + dble(Latt%listk(n,2))*Latt%b2_p  
      enddo
      
      ! Do timedisplaced k-resolved ===============================
      do nb = 1, nbins
         Phase (nb) = sgn(nb+n_skip)
         PhaseI(nb) = cmplx(sgn(nb+n_skip),0.d0,Kind(0.d0))
         do n = 1,Nunit
            do nt = 1,Lt_eff
               do no = 1,Norb
                  if (PartHole) then
                     bins(n,nt,nb) = bins(n,nt,nb) + ( bins_raw(n,nt,no,no,nb+n_skip) + bins_raw(n,Lt-nt+1,no,no,nb+n_skip) ) &
                                          & / cmplx(2.d0,0.d0,Kind(0.d0))
                  else
                     bins(n,nt,nb) = bins(n,nt,nb) + bins_raw(n,nt,no,no,nb+n_skip)
                  endif
                  if ( sqrt(Xk_p(1,n)**2 + Xk_p(2,n)**2) < 1.D-6 .and. N_Back == 1 ) then
                     bins(n,nt,nb) = bins(n,nt,nb) - Bins0_raw(no,nb+n_skip)*Bins0_raw(no,nb+n_skip) &
                                                     & *cmplx(dble(Nunit)/sgn(nb+n_skip),0.d0,kind(0.d0))
                  endif
               enddo
            enddo
         enddo
      enddo
      
      do n = 1,Nunit
         if (  Xk_p(1,n) >= -zero .and. XK_p(2,n) >= -zero ) then
            call COV(bins(n,:,:), phase, Xcov, Xmean, N_rebin )
            write(File_out,'(A,"_",F4.2,"_",F4.2)')  trim(name_obs), Xk_p(1,n), Xk_p(2,n)
            Open (Unit=10,File=File_out,status="unknown")
            Write(10,*) Lt_eff, nbins/N_rebin, real(lt-1,kind(0.d0))*dtau
            do nt = 1, LT_eff
               Write(10,"(F14.7,2x,F16.8,2x,F16.8)") &
                     & dble(nt-1)*dtau,  dble(Xmean(nt)), sqrt(abs(dble(Xcov(nt,nt))))
            enddo
            If (N_cov == 1) Then ! print covarariance
               Do nt = 1,LT_eff
                  Do nt1 = 1,LT_eff
                     Write(10,*) dble(Xcov(nt,nt1))
                  Enddo
               Enddo
            Endif
            close(10)
         endif
      enddo
      
      
      ! Do timedisplaced r=0 ===============================
      V_help = cmplx(0.d0,0.d0,kind(0.d0))
      do n = 1,Nunit
         do nb = 1,nbins
            do nt = 1,LT_eff
               V_help(nt,nb) = V_help(nt,nb) + bins(n,nt,nb)
            enddo
         enddo
      enddo
      V_help = V_help/dble(Nunit)
      call COV(V_help, phase, Xcov, Xmean, N_Rebin )
      write(File_out,'(A,"_R0")') trim(name_obs)
      Open (Unit=10,File=File_out,status="unknown")
      Write(10,*) LT_eff, nbins/N_rebin, real(lt-1,kind(0.d0))*dtau
      do nt = 1, LT_eff
         Write(10,"(F14.7,2x,F16.8,2x,F16.8)") &
               & dble(nt-1)*dtau,  dble(Xmean(nt)), sqrt(abs(dble(Xcov(nt,nt))))
      enddo
      If (N_cov == 1) Then ! Print  covariance
         Do nt = 1,LT_eff
            Do nt1 = 1,LT_eff
               Write(10,*) dble(Xcov(nt,nt1))
            Enddo
         Enddo
      Endif
      close(10)
      
      
      ! Do susceptibilities ===============================
      do nb = 1, nbins
         do n = 1,Nunit
            Z = cmplx(0.d0,0.d0,kind(0.d0))
            Do nt = 1,Lt_eff -1
               do no = 1,Norb
                  do no1 = 1,Norb
                     Z = Z + cmplx(0.5d0,0.d0,Kind(0.d0)) * ( bins_raw(n,nt,no,no1,nb+n_skip) + bins_raw(n,nt+1,no,no1,nb+n_skip) )
                     if ( sqrt(Xk_p(1,n)**2 + Xk_p(2,n)**2) < 1.D-6 .and. N_Back == 1 ) then
                        z = z - Bins0_raw(no,nb+n_skip)*Bins0_raw(no1,nb+n_skip) * cmplx(dble(Nunit)/sgn(nb+n_skip),0.d0,kind(0.d0))
                     endif
                  enddo
               enddo
            enddo
            if (PartHole) Z = Z*cmplx(2.d0,0.d0,Kind(0.d0))
            Bins_chi(N,nb) = Z 
         enddo
      enddo
      write(File_out,'(A,"_tauJK")') trim(name_obs)
      Open (Unit=33,File=File_out ,status="unknown")
      Do n = 1,Nunit
         call ERRCALCJ(Bins_chi(n,:), PhaseI, ZMean, ZERR, N_rebin )
         Zmean = Zmean*dtau
         Zerr = Zerr*dtau
         Write(33,"(F12.6,2x,F12.6,2x,F16.8,2x,F16.8,2x,F16.8,2x,F16.8)") &
               &   Xk_p(1,n), Xk_p(2,n), dble(ZMean), dble(ZERR), aimag(ZMean), aimag(ZERR)
      enddo
      Close(33)
      
      
      ! Deallocate space ===============================
      Deallocate ( bins, Phase,PhaseI, Xk_p, V_help)
      Deallocate ( Bins_chi, Xmean, Xcov )
      
   end Subroutine ana_tau

!==============================================================================
   
   subroutine Cov_eq(file)
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Analysis program for equal time observables.
!> 
!
!--------------------------------------------------------------------

      Implicit none
      Character (len=64), intent(in) :: file
      
      Real    (Kind=Kind(0.d0)), allocatable :: sgn(:)
      Complex (Kind=Kind(0.d0)), pointer :: Bins_raw(:,:,:,:,:), Bins0_raw(:,:)
      Type (Lattice)   :: Latt
      Real    (Kind=Kind(0.d0)) :: dtau
      
      call read_latt(file, sgn, bins_raw, bins0_raw, Latt, dtau)
      call ana_eq(file, sgn, bins_raw, bins0_raw, Latt)
   
   end subroutine Cov_eq

!==============================================================================
   
   Subroutine ana_eq(name, sgn, bins_raw, bins0_raw, Latt)
      Implicit none
      Character (len=64), intent(in) :: name
      Real    (Kind=Kind(0.d0)), allocatable, intent(in) :: sgn(:)
      Complex (Kind=Kind(0.d0)), pointer, intent(in) :: Bins_raw(:,:,:,:,:)
      Complex (Kind=Kind(0.d0)), pointer, intent(in) :: Bins0_raw(:,:)
      Type (Lattice), intent(in)    :: Latt
      
      Character (len=64)  :: File_out
      Integer :: N_skip, N_rebin, N_Cov, N_Back, N_auto
      Integer :: Nbins, Norb, Nunit
      Integer :: i, n, nb, no, no1, ierr
      Type  (Mat_C), allocatable :: Bins (:,:), Bins_R(:,:)
      Complex (Kind=Kind(0.d0)), allocatable :: Phase(:)
      Complex (Kind=Kind(0.d0)), allocatable :: V_help(:), V_help_R(:)
      Complex (Kind=Kind(0.d0)), allocatable :: Bins0(:,:)
      Real (Kind=Kind(0.d0)), allocatable :: Xk_p_s(:,:)
      Real (Kind=Kind(0.d0)), allocatable :: AutoCorr(:),En(:)
      Real    (Kind=Kind(0.d0)) :: Xk_p(2), Xr_p(2)
      Complex (Kind=Kind(0.d0)) :: Z, Xmean, Xerr, Xmean_r, Xerr_r
      Real (Kind=Kind(0.d0)) :: Xm,Xe

      NAMELIST /VAR_errors/   N_skip, N_rebin, N_Cov, N_Back, N_auto
         
      N_Back = 1
      N_auto = 0
      OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'unable to open <parameters>',ierr
         STOP
      END IF
      READ(5,NML=VAR_errors)
      CLOSE(5)
      
      Nbins = size(bins_raw,5)
      Norb  = size(bins_raw,3)
      Nunit = Latt%N
      
      Write(6,*) "# of bins: ", Nbins
      Nbins  = Nbins - n_skip
      Write(6,*) "Effective # of bins: ", Nbins/N_rebin
      N_auto=min(N_auto,Nbins/3)
      if(Nbins/N_rebin < 2) then
         stop "Effective # of bins smaller than 2. Analysis impossible!"
      endif

      ! Allocate  space
      Allocate ( bins(Nunit,Nbins), bins_r(Nunit,Nbins), Phase(Nbins), V_help(Nbins), V_help_R(Nbins), Bins0(Nbins,Norb))
      Do n = 1,Nunit
         do nb = 1,nbins
            Call Make_Mat(bins  (n,nb),Norb)
            Call Make_Mat(bins_r(n,nb),Norb)
            bins_r(n,nb)%el = cmplx(0.d0,0.d0,kind(0.d0))
            bins  (n,nb)%el = cmplx(0.d0,0.d0,kind(0.d0))
         Enddo
      Enddo 
         
      allocate( Xk_p_s(2, Latt%N) )
      do n = 1,Latt%N
            Xk_p_s(:,n) = dble(Latt%listk(n,1))*Latt%b1_p + dble(Latt%listk(n,2))*Latt%b2_p
      enddo
      
      Bins0 = cmplx(0.d0,0.d0,kind(0.d0))
      do nb = 1, nbins + n_skip
         if (nb > n_skip ) then
            Phase(nb-n_skip) = cmplx(sgn(nb),0.d0,kind(0.d0))
            Do no = 1,Norb
               if (N_Back == 1 ) Bins0(nb-n_skip,no) = Bins0_raw(no,nb)
            enddo
            do n = 1,Nunit
               do no = 1,norb
                  do no1 = 1,Norb
                     bins(n,nb-n_skip)%el(no,no1) = Bins_raw(n,1,no,no1,nb)
                  enddo
               enddo
               Xk_p(:) = Xk_p_s(:,n)
               if ( sqrt(Xk_p(1)**2 + Xk_p(2)**2) < 1.D-6 ) then
                  do no = 1,norb
                     do no1 = 1,Norb
                        bins(n,nb-n_skip)%el(no,no1)  =  bins(n,nb-n_skip)%el(no,no1) -  &
                              &        cmplx(dble(Latt%N),0.d0,kind(0.d0))*Bins0(nb-n_skip,no)*Bins0(nb-n_skip,no1) &
                              &        /Phase(nb-n_skip)
                     enddo
                  enddo
               endif
            enddo
         endif
      enddo
      N_auto=min(N_auto,Nbins/3)
      

      Call Fourier_K_to_R(bins,bins_r,Latt)

      ! Setup symmetries for square lattice. 
#ifdef test
      do n = 1,Nunit
         n1 = n
         Write(6,*) Xk_p(1,n1), Xk_p(2,n1)
         do m = 1,4
            n1 = Rot90(n1, Xk_p, Nunit)
            Write(6,*) n1, Xk_p(1,n1), Xk_p(2,n1)
         enddo
         Write(6,*)
      enddo
#endif
      write(File_out,'(A,A)') trim(name), "JK"
      Open (Unit=33,File=File_out ,status="unknown")
      write(File_out,'(A,A)') trim(name), "JR"
      Open (Unit=34,File=File_out ,status="unknown")
      Do n = 1,Nunit
         Xk_p = dble(Latt%listk(n,1))*Latt%b1_p + dble(Latt%listk(n,2))*Latt%b2_p 
         Xr_p = dble(Latt%list (n,1))*Latt%a1_p + dble(Latt%list (n,2))*Latt%a2_p 
         Write(33,"(F12.6,2x,F12.6)")  Xk_p(1), Xk_p(2)
         Write(34,"(F12.6,2x,F12.6)")  Xr_p(1), Xr_p(2)
         Do no = 1,Norb
            do no1 = 1,Norb
               do nb = 1,Nbins
                  V_help(nb) = bins  (n,nb)%el(no,no1)
               enddo
               call ERRCALCJ( V_help, Phase,XMean, XERR, N_rebin ) 
               Write(33,"(I3,2x,I3,2x,F16.8,2x,F16.8,2x,F16.8,2x,F16.8)") &
                     &  no,no1, dble(XMean), dble(XERR), aimag(XMean), aimag(XERR)
               do nb = 1,Nbins
                  V_help(nb) = bins_r(n,nb)%el(no,no1)
               enddo
               call ERRCALCJ( V_help,Phase, XMean_r, XERR_r, N_rebin ) 
               Write(34,"(I3,2x,I3,2x,F16.8,2x,F16.8,2x,F16.8,2x,F16.8)") &
                     &  no,no1, dble(XMean_r), dble(XERR_r), aimag(XMean_r), aimag(XERR_r)
            enddo
         enddo
      enddo

      Close(33)
      Close(34)
      
      if(N_auto>0) then
      ALLOCATE(AutoCorr(N_auto))
      ALLOCATE(EN(Nbins))
      Do n = 1,Nunit
         Xk_p = dble(Latt%listk(n,1))*Latt%b1_p + dble(Latt%listk(n,2))*Latt%b2_p
         if (Xk_p(1) >= -1.d-8 .and. XK_p(2) >= -1.d-8) then
            write(File_out,'(A,"_Auto_Tr_",F4.2,"_",F4.2)') trim(name), Xk_p(1), Xk_p(2)
            OPEN (UNIT=21, FILE=File_out, STATUS='unknown')
            WRITE(21,*)
            do nb = 1,Nbins
               Z=0
               do no = 1,Norb
                  Z = Z+bins  (n,nb)%el(no,no)
               enddo
               En(nb)=dble(Z)
            enddo
            Call AUTO_COR(En,AutoCorr)
            do i = 1,N_auto
               CALL ERRCALCJ(En,XM, XE,i)
               write(21,*) i, AutoCorr(i), Xe
            enddo
            CLOSE(21)
         endif
      enddo
      DEALLOCATE(AutoCorr)
      endif
      
   end Subroutine ana_eq


   subroutine Cov_vec(file)
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Analysis program for scalar observables
!
!--------------------------------------------------------------------
        
      Implicit none
      Character (len=64), intent(in) :: file

      Real    (Kind=Kind(0.d0)), allocatable :: sgn_raw(:)
      Complex (Kind=Kind(0.d0)), pointer     :: Bins_raw(:,:)
        
      call read_vec(file, sgn_raw, bins_raw)
      call ana_vec(file, sgn_raw, bins_raw)
      
   END subroutine Cov_vec

!==============================================================================

   subroutine ana_vec(name, sgn_raw, bins_raw)
      Implicit none
      Character (len=64), intent(in) :: name
      Real    (Kind=Kind(0.d0)), allocatable, intent(inout) :: sgn_raw(:)
      Complex (Kind=Kind(0.d0)), pointer,     intent(inout) :: bins_raw(:,:)

      REAL    (Kind=Kind(0.d0)), DIMENSION(:),   ALLOCATABLE :: EN, sgn
      REAL    (Kind=Kind(0.d0)) :: XM, XERR, X

      Complex (Kind=Kind(0.d0)), Allocatable  :: Bins(:,:)
      REAL    (Kind=Kind(0.d0)), Allocatable  :: AutoCorr(:)
      Integer :: Nobs 
      Integer :: Nbins, Nbins_eff, I, IOBS, N_Back
      Integer :: N,N1

      Integer :: N_skip, N_rebin, N_Cov, ierr, N_auto
      Character (len=64) :: File_out
      NAMELIST /VAR_errors/   N_skip, N_rebin, N_Cov, N_Back, N_auto
      
      !New Stuff for Autocorrelation
      Integer :: tau, tau_max, Nbins_eff2
      REAL(Kind=Kind(0.d0)), DIMENSION(:)  , ALLOCATABLE :: vec, vec_err

      N_auto=0
      OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) 'unable to open <parameters>',ierr
         STOP
      END IF
      READ(5,NML=VAR_errors)
      CLOSE(5)
        
      Nobs  = size(bins_raw, 1)
      Nbins = size(bins_raw, 2)

      Write(6,*) "# of bins: ", Nbins
      Nbins_eff  = Nbins - n_skip
      Write(6,*) "Effective # of bins: ", Nbins_eff/N_rebin
      N_auto=min(N_auto,Nbins_eff/3)
      if(Nbins_eff/N_rebin < 2) then
         stop "Effective # of bins smaller than 2. Analysis impossible!"
      endif
        
      ! Allocate  space
      Allocate ( Bins(Nobs,Nbins_eff), sgn(Nbins_eff) )
      
      do i =1,Nbins_eff
         Bins(:,i) = Bins_raw(:,i+n_skip)
         sgn(i) = sgn_raw(i+n_skip)
      enddo
         
      write(File_out,'(A,A)') trim(name), "J"
      OPEN (UNIT=21, FILE=File_out, STATUS='unknown')
      WRITE(21,*) 'Effective number of bins, and bins: ', Nbins_eff/N_rebin, Nbins
      ALLOCATE (EN(Nbins_eff), vec(NOBS), vec_err(NOBS))
      DO IOBS = 1,NOBS
         EN(:) = Real(Bins(IOBS,:), kind(0.d0))
         CALL ERRCALCJ(EN,sgn,XM,XERR,N_Rebin)
         vec    (IOBS) = XM
         vec_err(IOBS) = XERR
         WRITE(21,*)
         WRITE(21,2001) IOBS, XM,  XERR
      ENDDO
      CALL ERRCALCJ(sgn, XM,XERR,N_Rebin)
      WRITE(21,*)
      WRITE(21,2001) NOBS+1, XM,  XERR
      CLOSE(21)
2001    FORMAT('OBS : ', I4,4x,F12.6,2X, F12.6)
!2001  FORMAT('OBS : ', I4,4x,ES12.5,2X, ES12.5)
        
      if(N_auto>0) then
         ALLOCATE(AutoCorr(N_auto))
         DO IOBS = 1,NOBS
            write(File_out,'(A,A,I1.1)') trim(name), '_Auto_', iobs
            write(*,*) File_out
            OPEN (UNIT=21, FILE=File_out, STATUS='unknown')
            WRITE(21,*)
            EN(:) = Real(Bins(IOBS,:), kind(0.d0))
            Call AUTO_COR(EN,AutoCorr)
            do i = 1,N_auto
               CALL ERRCALCJ(EN,XM,XERR,i)
               write(21,*) i, AutoCorr(i), Xerr
            enddo
            CLOSE(21)
         ENDDO
         DEALLOCATE(AutoCorr)
      endif

      DEALLOCATE (EN,vec,vec_err,sgn_raw,sgn,Bins_raw,Bins)
      
   END subroutine ana_vec
end module ana_mod
