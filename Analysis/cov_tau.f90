!  Copyright (C) 2016 The ALF project
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

       Program Cov_tau

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Analysis of imaginary time displaced  correlation functions.
!
!--------------------------------------------------------------------         
         Use Errors
         Use MyMats
         Use Matrix

         Implicit none


         Integer :: Nunit, Norb, N_auto
         Integer :: no, no1, n, nbins, n_skip, nb, NT, NT1, Lt, N_rebin, N_cov, ierr, N_Back
         real    (Kind=Kind(0.d0)):: X, Y,  dtau, X_diag
         Complex (Kind=Kind(0.d0)), allocatable :: Xmean(:), Xcov(:,:)
         Complex (Kind=Kind(0.d0)) :: Z, Z_diag
         Real    (Kind=Kind(0.d0)) :: Zero=1.D-8
         Real    (Kind=Kind(0.d0)), allocatable :: Phase(:)
         Complex (Kind=Kind(0.d0)), allocatable :: Bins(:,:,:)
         Complex (Kind=Kind(0.d0)), allocatable :: Bins0(:,:)
         Complex (Kind=Kind(0.d0)), allocatable :: V_help(:,:)
         Real    (Kind=Kind(0.d0)), allocatable :: Xk_p(:,:)
         Character (len=64) :: File_out

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

         ! Determine the number of bins. 
         Open ( Unit=10, File="intau", status="unknown" ) 
         nbins = 0
         do
            Read(10,*,End=10) X,Norb,Nunit, LT, dtau
            Do no = 1,Norb
               Read(10,*) Z
            enddo
            do n = 1,Nunit
               Read(10,*) X,Y
               do nt = 1,LT
                  do no = 1,Norb
                     do no1 = 1,Norb
                        read(10,*) Z
                     enddo
                  enddo
               enddo
            enddo
            !Write(6,*) nbins
            nbins = nbins + 1
         enddo
10       continue
         Close(10) 
         Write(6,*) "# of bins: ", Nbins
         nbins  = Nbins - n_skip
         Write(6,*) "Effective # of bins: ", Nbins


         ! Allocate  space
         Allocate ( bins(Nunit,Lt,Nbins), Phase(Nbins), Xk_p(2,Nunit), V_help(lt,Nbins), bins0(Nbins,Norb))
         
         Allocate (Xmean(Lt), Xcov(Lt,Lt))
         bins  = 0.d0
         bins0 = cmplx(0.d0,0.d0,Kind(0.d0))
         Open ( Unit=10, File="intau", status="unknown" ) 
         do nb = 1, nbins + n_skip
            if (nb > n_skip ) then
               Read(10,*,End=10) Phase(nb-n_skip),no,no1,n, X
               Z_diag = cmplx(0.d0,0.d0,kind(0.d0))
               Do no = 1,Norb
                  Read(10,*)   Z 
                  If ( N_back == 1 )   bins0(nb-n_skip,no) = Z
                  Z_diag =  Z_diag + bins0(nb-n_skip,no)*bins0(nb-n_skip,no)
               Enddo
               do n = 1,Nunit
                  Read(10,*) Xk_p(1,n), Xk_p(2,n)
                  do nt = 1,Lt
                     do no = 1,norb
                        do no1 = 1,Norb
                           read(10,*) Z
                           if (no == no1) bins(n,nt,nb-n_skip) = bins(n,nt,nb-n_skip) +  Z
                        enddo
                     enddo
                  enddo
                  if ( sqrt(Xk_p(1,n)**2 + Xk_p(2,n)**2) < 1.D-6 ) then
                     Do nt = 1,Lt
                        bins(n,nt,nb-n_skip) =   bins(n,nt,nb-n_skip) - Z_diag*cmplx(dble(Nunit),0.d0,kind(0.d0))
                     enddo
                  endif
               enddo
            else
               Read(10,*,End=10) X,no,no1,n,Y
               do no = 1,Norb
                  Read(10,*) Z
               enddo
               do n = 1,Nunit
                  Read(10,*) X,Y
                  do nt = 1,LT
                     do no = 1,Norb
                        do no1 = 1,Norb
                           read(10,*) Z
                        enddo
                     enddo
                  enddo
               enddo
            endif
         enddo
         close(10)

         !do n = 1,Nbins
         !   Write(6,*) Phase(n)
         !Enddo
         do n = 1,Nunit
            if (  Xk_p(1,n) >= -zero .and. XK_p(2,n) >= -zero ) then
               call COV(bins(n,:,:), phase, Xcov, Xmean, N_rebin )
               write(File_out,'("g_",F4.2,"_",F4.2)')  Xk_p(1,n), Xk_p(2,n)
               Open (Unit=10,File=File_out,status="unknown")
               do nt = 1, LT
                  Write(10,"(F14.7,2x,F16.8,2x,F16.8)") &
                       & dble(nt-1)*dtau,  dble(Xmean(nt)), sqrt(abs(dble(Xcov(nt,nt))))
               enddo
               If (N_cov == 1) Then ! print covarariance
                  Do nt = 1,LT
                     Do nt1 = 1,LT
                        Write(10,*) dble(Xcov(nt,nt1))
                     Enddo
                  Enddo
               Endif
               close(10)
            endif
         enddo

         V_help = cmplx(0.d0,0.d0,kind(0.d0))
         do n = 1,Nunit
            do nb = 1,nbins
               do nt = 1,LT
                  V_help(nt,nb) = V_help(nt,nb) +  bins(n,nt,nb)
               enddo
            enddo
         enddo
         V_help = V_help/dble(Nunit)
         call COV(V_help, phase, Xcov, Xmean, N_Rebin )
         write(File_out,'("g_R0")') 
         Open (Unit=10,File=File_out,status="unknown")
         do nt = 1, LT
            Write(10,"(F14.7,2x,F16.8,2x,F16.8)") &
                 & dble(nt-1)*dtau,  dble(Xmean(nt)), sqrt(abs(dble(Xcov(nt,nt))))
         enddo
         If (N_cov == 1) Then ! Print  covariance
            Do nt = 1,LT
               Do nt1 = 1,LT
                  Write(10,*) dble(Xcov(nt,nt1))
               Enddo
            Enddo
         Endif
         close(10)
         

       end Program Cov_tau


!!$         Interface
!!$            Integer function Rot90(n, Xk_p, Nunit)
!!$              Implicit none
!!$              Integer, INTENT(IN)       :: Nunit,n
!!$              Real (Kind=Kind(0.d0)), INTENT(IN) :: Xk_p(2,Nunit)
!!$            end function Rot90
!!$         end Interface
!!$
!!$       Integer function Rot90(n, Xk_p, Nunit)
!!$         
!!$         Implicit none
!!$         Integer, INTENT(IN)       :: Nunit,n
!!$         Real (Kind=Kind(0.d0)), INTENT(IN) :: Xk_p(2,Nunit)
!!$         
!!$         !Local
!!$         real (Kind=Kind(0.d0)) :: X1_p(2), Zero, pi, X
!!$         Integer :: m
!!$         
!!$         Zero = 1.D-4
!!$         pi = acos(-1.d0)
!!$         X1_p(1)  =  Xk_p(2,n)   
!!$         X1_p(2)  = -Xk_p(1,n)   
!!$         if (X1_p(1) < -pi + Zero )  X1_p(1) = X1_p(1) + 2.0*pi
!!$         if (X1_p(2) < -pi + Zero )  X1_p(2) = X1_p(2) + 2.0*pi
!!$         
!!$         Rot90 = 0
!!$         Do m = 1,Nunit
!!$            X = sqrt( (X1_p(1) -Xk_p(1,m))**2 +  (X1_p(2) -Xk_p(2,m))**2 )
!!$            If ( X < Zero) then
!!$               Rot90 = m
!!$               exit
!!$            endif
!!$         Enddo
!!$         
!!$       end function Rot90
