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


   Program ana
      use ana_mod
      implicit none
      Integer                         :: i, n, nargs
      Character (len=64)              :: name
      Character (len=64), allocatable :: names(:)
      Logical                         :: PartHole
      
      nargs = COMMAND_ARGUMENT_COUNT()
      allocate( names(nargs) )
      do i = 1, nargs
         CALL GET_COMMAND_ARGUMENT(i, name)
         names(i) = name
      enddo
      
      do n=1, size(names)
         name = names(n)
         i = len(trim(name)) -4
         if ( name(i:) == '_scal' ) then
            print *, ''
            print '(A,A)', "analyzing ", name
            call Cov_vec(name)
         endif
      enddo
      
      do n=1, size(names)
         name = names(n)
         i = len(trim(name)) -2
         if ( name(i:) == '_eq' ) then
            print *, ''
            print '(A,A)', "analyzing ", name
            call Cov_eq(name)
         endif
      enddo
      
      do n=1, size(names)
         name = names(n)
         i = len(trim(name)) -3
         if ( name(i:) == '_tau' ) then
            if (name == 'Green_tau') then
               PartHole = .false.
            else
               PartHole = .true.
            endif
            print *, ''
            print '(A,A)', "analyzing ", name
            print *, PartHole
            call Cov_tau(name, PartHole)
         endif
      enddo
      
   end Program ana
       
