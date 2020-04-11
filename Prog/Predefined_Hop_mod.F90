!  Copyright (C) 2016 - 2020 The ALF project
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
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
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

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>
!> @brief 
!> This module provides a set of predefined hoppings as well as a 
!> general framework to specify the hopping matrix for translation invariant      
!> multi-orbital systems.
!> 
!--------------------------------------------------------------------

    Module Predefined_Hoppings
      
      Use Lattices_v3
      Use Operator_mod
      Use WaveFunction_mod
      Use MyMats
      Implicit none

      
      Type Hopping_Matrix_type
         Integer                   :: N_bonds
         Complex (Kind=Kind(0.d0)), pointer :: T    (:)    !  This does not include  local terms
         Complex (Kind=Kind(0.d0)), pointer :: T_loc(:)    !  This is just for the local matrix elements such as chemical potential-
         Integer                  , pointer :: list(:,:)
         ! T(N_b=1..N_bonds)
         ! List(N_b,1) = no_1
         ! List(N_b,2) = no_2
         ! List(N_b,3) = n_1
         ! List(N_b,4) = n_2
         ! H_[(i,no_1),(i + n_1 a_1 + n_2 a_2,no_2)] = T(N_b) 
         Integer                   :: N_Phi   
         Real    (Kind=Kind(0.d0)) :: Phi_X, Phi_Y
         Logical                   :: Bulk
         ! N_Phi         = #  of flux quanta  piercieng the lattice
         ! Phi_X, Phi_Y  =  Twist
         ! Bulk          =  Twist as boundary condtion (Bulk=.F.) 
         !               =  Twist as Vector potential  (Bulk=.T.)
      end type Hopping_Matrix_Type

      Type (Hopping_Matrix_type), allocatable, private ::  Hopping_Matrix(:)
      Integer, private              ::  N_Fam
      Integer, allocatable, private ::  L_Fam(:),  List_Fam(:,:,:), Multiplicity(:)
      Real (Kind=Kind(0.d0)), allocatable,  private ::  Prop_Fam(:)
      
      
    contains

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>
!> @brief 
!> Default hopping for the square lattice.  Ham_T is the nearest neighbour hopping and Ham_Chem the chemical potention.
!> 
!     
!--------------------------------------------------------------------      
      Subroutine Set_Default_hopping_parameters_square(Ham_T, Ham_Chem, Phi_X, Phi_Y, Bulk,  N_Phi, N_FL, &
           &                                           List, Invlist, Latt, Latt_unit, Checkerboard )
 
        Implicit none
        
        Real (Kind=Kind(0.d0)), Intent(IN)    :: Ham_T, Ham_Chem, Phi_x, Phi_y
        Integer, Intent(IN)                   :: N_Phi, N_FL
        Logical, Intent(IN)                   :: Bulk
        Integer, Intent(IN), Dimension(:,:)   :: List, Invlist
        Type(Lattice),  Intent(in)            :: Latt
        Type(Unit_cell),Intent(in)            :: Latt_unit
        Logical,  Intent (in)                 :: Checkerboard


        ! Local
        Integer :: nf,N_Bonds, nc, I, I1
        Real (Kind=Kind(0.d0)) :: Zero = 1.0E-8

        Allocate( Hopping_Matrix(N_FL) )
        do nf = 1,N_FL
           Hopping_Matrix(nf)%N_bonds = 0
           if ( abs(Ham_T) > Zero)  then
              Hopping_Matrix(nf)%N_bonds = 2
              Allocate (Hopping_Matrix(nf)%List(Hopping_Matrix(nf)%N_bonds,4), &
                   &    Hopping_Matrix(nf)%T(Hopping_Matrix(nf)%N_bonds) )
              nc = 0
              nc = nc + 1
              Hopping_Matrix(nf)%T(nc)    = cmplx(-Ham_T,0.d0,kind(0.d0))
              Hopping_Matrix(nf)%List(nc,1) = 1
              Hopping_Matrix(nf)%List(nc,2) = 1
              Hopping_Matrix(nf)%List(nc,3) = 0
              Hopping_Matrix(nf)%List(nc,4) = 1
              
              nc = nc + 1
              Hopping_Matrix(nf)%T(nc)    = cmplx(-Ham_T,0.d0,kind(0.d0))
              Hopping_Matrix(nf)%List(nc,1) = 1
              Hopping_Matrix(nf)%List(nc,2) = 1
              Hopping_Matrix(nf)%List(nc,3) = 1
              Hopping_Matrix(nf)%List(nc,4) = 0
           Endif
           Allocate ( Hopping_Matrix(nf)%T_Loc(Latt_Unit%Norb) )
           do nc = 1,Latt_Unit%Norb
              Hopping_Matrix(nf)%T_Loc(nc)  = cmplx(-Ham_Chem,0.d0,kind(0.d0))
           enddo
           Hopping_Matrix(nf)%N_Phi =  N_Phi
           Hopping_Matrix(nf)%Phi_X =  Phi_X
           Hopping_Matrix(nf)%Phi_Y =  Phi_Y
           Hopping_Matrix(nf)%Bulk =   Bulk
        enddo

        If  (Checkerboard)  then
           N_Fam  = 4
           Allocate (L_FAM(N_FAM),  Prop_Fam(N_FAM), Multiplicity(Latt_unit%Norb) )
           L_FAM  = Latt%N/2
           Prop_Fam = 1.d0
           Allocate (List_Fam(N_Fam,L_FAM(1),2))
           Multiplicity = 4
           L_FAM  = 0
           do I = 1,Latt%N
              if ( mod(Latt%List(I,1) + Latt%List(I,2),2) == 0 ) then
                 Nf = 1
                 L_FAM(Nf) = L_FAM(Nf) + 1
                 List_Fam(Nf,L_FAM(Nf),1) = I ! Unit cell
                 List_Fam(Nf,L_FAM(Nf),2) = 1 ! The bond (See above)
                 Nf = 2
                 L_FAM(Nf) = L_FAM(Nf) + 1
                 List_Fam(Nf,L_FAM(Nf),1) = I
                 List_Fam(Nf,L_FAM(Nf),2) = 2 
              else
                 Nf = 3
                 L_FAM(Nf) = L_FAM(Nf) + 1
                 List_Fam(Nf,L_FAM(Nf),1) = I
                 List_Fam(Nf,L_FAM(Nf),2) = 1  
                 Nf = 4
                 L_FAM(Nf) = L_FAM(Nf) + 1
                 List_Fam(Nf,L_FAM(Nf),1) = I
                 List_Fam(Nf,L_FAM(Nf),2) = 2  
              endif
           enddo
        endif
        
      end Subroutine Set_Default_hopping_parameters_square

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>
!> @brief 
!> Default hopping for n-leg-ladder.  Ham_T is the nearest neighbour hopping along the chain,  Ham_T_perp  the 
!> interrung hopping and H_chem the chemical potential.
!> 
!     
!--------------------------------------------------------------------      
      Subroutine Set_Default_hopping_parameters_N_Leg_Ladder(Ham_T, Ham_T_perp,  Ham_Chem, Phi_X, Phi_Y, Bulk,  N_Phi, N_FL, &
           &                                           List, Invlist, Latt, Latt_unit, Checkerboard )
 
        Implicit none
        
        Real (Kind=Kind(0.d0)), Intent(IN)    :: Ham_T, Ham_T_perp, Ham_Chem, Phi_x, Phi_y
        Integer, Intent(IN)                   :: N_Phi, N_FL
        Logical, Intent(IN)                   :: Bulk
        Integer, Intent(IN), Dimension(:,:)   :: List, Invlist
        Type(Lattice),  Intent(in)            :: Latt
        Type(Unit_cell),Intent(in)            :: Latt_unit
        Logical,  Intent (in)                 :: Checkerboard


        ! Local
        Integer :: nf,N_Bonds, nc, I, I1, n, no
        Real (Kind=Kind(0.d0)) :: Zero = 1.0E-8

        Allocate( Hopping_Matrix(N_FL) )
        do nf = 1,N_FL
           Hopping_Matrix(nf)%N_bonds = Latt_unit%Norb +  (Latt_unit%Norb - 1 )
           Allocate (Hopping_Matrix(nf)%List(Hopping_Matrix(nf)%N_bonds,4), &
                &    Hopping_Matrix(nf)%T(Hopping_Matrix(nf)%N_bonds) )
           nc = 0
           do n = 1,Latt_unit%Norb
              nc = nc + 1
              Hopping_Matrix(nf)%T(nc)    = cmplx(-Ham_T,0.d0,kind(0.d0))
              Hopping_Matrix(nf)%List(nc,1) = n
              Hopping_Matrix(nf)%List(nc,2) = n
              Hopping_Matrix(nf)%List(nc,3) = 1
              Hopping_Matrix(nf)%List(nc,4) = 0
           enddo
           
           do n = 1,Latt_unit%Norb -1 
              nc = nc + 1
              Hopping_Matrix(nf)%T(nc)    = cmplx(-Ham_T_perp,0.d0,kind(0.d0))
              Hopping_Matrix(nf)%List(nc,1) = n
              Hopping_Matrix(nf)%List(nc,2) = n + 1
              Hopping_Matrix(nf)%List(nc,3) = 0
              Hopping_Matrix(nf)%List(nc,4) = 0
           enddo
           
           Allocate ( Hopping_Matrix(nf)%T_Loc(Latt_Unit%Norb) )
           do nc = 1,Latt_Unit%Norb
              Hopping_Matrix(nf)%T_Loc(nc)  = cmplx(-Ham_Chem,0.d0,kind(0.d0))
           enddo
           Hopping_Matrix(nf)%N_Phi =  N_Phi
           Hopping_Matrix(nf)%Phi_X =  Phi_X
           Hopping_Matrix(nf)%Phi_Y =  Phi_Y
           Hopping_Matrix(nf)%Bulk =   Bulk
        enddo

        !Write(6,*) Latt_unit%Norb
        If (Checkerboard)  then
           Allocate ( Multiplicity(Latt_unit%Norb) )
           If     ( Latt_Unit%Norb  == 1 ) then
              Multiplicity = 2
              N_Fam        = 2
           elseif ( Latt_Unit%Norb  == 2 ) then
              Multiplicity = 3
              N_Fam        = 3
           else
              Multiplicity                 = 4
              Multiplicity(1)              = 3
              Multiplicity(Latt_unit%Norb) = 3
              N_Fam        = 4
           endif
           Allocate ( L_FAM(N_FAM),  Prop_Fam(N_FAM) )
           L_Fam    = Latt%N*Latt_unit%Norb/2
           Prop_Fam = 1.d0
           Allocate ( List_Fam(N_Fam,L_FAM(1),2) )
           

           L_FAM  = 0
           do I = 1,Latt%N
              if ( mod(Latt%List(I,1),2) == 0 ) then
                 Nf = 1
                 do no = 1,Latt_unit%Norb
                    L_FAM(Nf) = L_FAM(Nf) + 1
                    List_Fam(Nf,L_FAM(Nf),1) = I ! Unit cell
                    List_Fam(Nf,L_FAM(Nf),2) = no ! The bond (See above)
                 enddo
              else 
                 Nf = 2
                 do no = 1,Latt_unit%Norb
                    L_FAM(Nf) = L_FAM(Nf) + 1
                    List_Fam(Nf,L_FAM(Nf),1) = I
                    List_Fam(Nf,L_FAM(Nf),2) = no
                 enddo
              endif
           enddo
           do no = 1,Latt_unit%Norb - 1
              if (mod(no,2) == 1 ) then
                 Nf = 3
                 Write(6,*)  NF, no + Latt_unit%Norb
                 do I = 1,Latt%N
                    L_FAM(Nf) = L_FAM(Nf) + 1
                    List_Fam(Nf,L_FAM(Nf),1) = I
                    List_Fam(Nf,L_FAM(Nf),2) = no + Latt_unit%Norb
                 enddo
              else
                 Nf = 4
                 Write(6,*)  NF, no + Latt_unit%Norb
                 do I = 1,Latt%N
                    L_FAM(Nf) = L_FAM(Nf) + 1
                    List_Fam(Nf,L_FAM(Nf),1) = I
                    List_Fam(Nf,L_FAM(Nf),2) = no + Latt_unit%Norb
                 enddo
              endif
           enddo
        endif
        
      end Subroutine Set_Default_hopping_parameters_N_Leg_Ladder

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>
!> @brief 
!> Default hopping for Honeycomb lattice.  Ham_T is the nearest neighbour hopping along the chain, Lambda is  the 
!> Kane-Mele term and H_chem the chemical potential.  **Note**  The Kane-Mele term is not yet implemented.
!> 
!     
!--------------------------------------------------------------------      
      Subroutine Set_Default_hopping_parameters_honeycomb(Ham_T, Ham_Lambda, Ham_Chem, Phi_X, Phi_Y, Bulk,  N_Phi, N_FL,&
                                                        & List, Invlist, Latt, Latt_unit, Checkerboard )
        
        Implicit none
        
        Real (Kind=Kind(0.d0)), Intent(IN)    :: Ham_T, Ham_Lambda, Ham_Chem, Phi_x, Phi_y
        Integer, Intent(IN)                   :: N_Phi, N_FL
        Logical, Intent(IN)                   :: Bulk
        Integer, Intent(IN), Dimension(:,:)   :: List, Invlist
        Type(Lattice),  Intent(in)            :: Latt
        Type(Unit_cell),Intent(in)            :: Latt_unit
        Logical,  Intent (in)                 :: Checkerboard


        ! Local
        Integer :: nf,N_Bonds, nc, I, I1, n, no
        Real (Kind=Kind(0.d0)) :: Zero = 1.0E-8

        If (abs(Ham_Lambda) > 0 ) then
           Write(6,*)  'Kane Mele term is not yet implemented'
           Stop
        endif
        Allocate( Hopping_Matrix(N_FL) )
        do nf = 1,N_FL
           Hopping_Matrix(nf)%N_bonds =  3
           Allocate (Hopping_Matrix(nf)%List(Hopping_Matrix(nf)%N_bonds,4), &
                &    Hopping_Matrix(nf)%T(Hopping_Matrix(nf)%N_bonds) )
           nc = 0
           nc = nc + 1
           Hopping_Matrix(nf)%T(nc)    = cmplx(-Ham_T,0.d0,kind(0.d0))
           Hopping_Matrix(nf)%List(nc,1) =  1
           Hopping_Matrix(nf)%List(nc,2) =  2
           Hopping_Matrix(nf)%List(nc,3) =  0
           Hopping_Matrix(nf)%List(nc,4) =  0
           
           nc = nc + 1
           Hopping_Matrix(nf)%T(nc)    = cmplx(-Ham_T,0.d0,kind(0.d0))
           Hopping_Matrix(nf)%List(nc,1) =  2
           Hopping_Matrix(nf)%List(nc,2) =  1 
           Hopping_Matrix(nf)%List(nc,3) =  0
           Hopping_Matrix(nf)%List(nc,4) =  1

           nc = nc + 1
           Hopping_Matrix(nf)%T(nc)    = cmplx(-Ham_T,0.d0,kind(0.d0))
           Hopping_Matrix(nf)%List(nc,1) =  1
           Hopping_Matrix(nf)%List(nc,2) =  2 
           Hopping_Matrix(nf)%List(nc,3) =  1
           Hopping_Matrix(nf)%List(nc,4) = -1

           Allocate ( Hopping_Matrix(nf)%T_Loc(Latt_Unit%Norb) )
           do nc = 1,Latt_Unit%Norb
              Hopping_Matrix(nf)%T_Loc(nc)  = cmplx(-Ham_Chem,0.d0,kind(0.d0))
           enddo
           Hopping_Matrix(nf)%N_Phi =  N_Phi
           Hopping_Matrix(nf)%Phi_X =  Phi_X
           Hopping_Matrix(nf)%Phi_Y =  Phi_Y
           Hopping_Matrix(nf)%Bulk =   Bulk
        enddo
        
        If  (Checkerboard)  then
           N_Fam  = 3
           Allocate (L_FAM(N_FAM),  Prop_Fam(N_FAM), Multiplicity(Latt_unit%Norb) )
           L_FAM  = Latt%N
           Prop_Fam = 1.d0
           Allocate (List_Fam(N_Fam,L_FAM(1),2))
           Multiplicity = 3
           do I = 1,Latt%N
              Do  nf = 1,N_Fam
                 List_Fam(nf,I,1) = I  ! Unit cell
                 List_Fam(nf,I,2) = nf ! The bond (See above)
              Enddo
           enddo
        endif
        
      end Subroutine Set_Default_hopping_parameters_honeycomb

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>
!> @brief 
!> Default hopping for a bilayer square. Ham_T1, Ham_T2, are the nearest neighbour hopping on the first and second layer and
!> Ham_T_perp   is the interlayer  hopping.
!> 
!     
!--------------------------------------------------------------------      
      Subroutine Set_Default_hopping_parameters_Bilayer_square(Ham_T1,Ham_T2,Ham_Tperp, Ham_Chem, Phi_X, Phi_Y, Bulk,  N_Phi, N_FL,&
           &                                                   List, Invlist, Latt, Latt_unit, Checkerboard )
 
        Implicit none
        
        Real (Kind=Kind(0.d0)), Intent(IN)    :: Ham_T1, Ham_T2, Ham_Tperp, Ham_Chem, Phi_x, Phi_y
        Integer, Intent(IN)                   :: N_Phi, N_FL
        Logical, Intent(IN)                   :: Bulk
        Integer, Intent(IN), Dimension(:,:)   :: List, Invlist
        Type(Lattice),  Intent(in)            :: Latt
        Type(Unit_cell),Intent(in)            :: Latt_unit
        Logical,  Intent (in)                 :: Checkerboard


        ! Local
        Integer :: nf,N_Bonds, nc, I, I1, No_Shift, n, nb
        Real (Kind=Kind(0.d0)) :: Zero = 1.0E-8
        Logical :: Test=.false.


        If (abs(Ham_T1) < Zero ) Then
           Write(6,*) 'At least Ham_T1 has to be bigger than zero'
           Stop
        endif

        
        Allocate( Hopping_Matrix(N_FL) )
        do nf = 1,N_FL
           N_bonds = 0
           if (abs(Ham_T1)    > Zero )  N_bonds = N_bonds + 2
           if (abs(Ham_Tperp) > Zero )  N_bonds = N_bonds + 1
           if (abs(Ham_T2)    > Zero )  N_bonds = N_bonds + 2
           Hopping_Matrix(nf)%N_bonds = N_bonds
           Allocate (Hopping_Matrix(nf)%List(Hopping_Matrix(nf)%N_bonds,4), &
                &    Hopping_Matrix(nf)%T(Hopping_Matrix(nf)%N_bonds) )
           nc = 0
           If (abs(Ham_T1) > Zero ) Then
              nc = nc + 1
              Hopping_Matrix(nf)%T(nc)    = cmplx(-Ham_T1,0.d0,kind(0.d0))
              Hopping_Matrix(nf)%List(nc,1) = 1
              Hopping_Matrix(nf)%List(nc,2) = 1
              Hopping_Matrix(nf)%List(nc,3) = 0
              Hopping_Matrix(nf)%List(nc,4) = 1
              
              nc = nc + 1
              Hopping_Matrix(nf)%T(nc)    = cmplx(-Ham_T1,0.d0,kind(0.d0))
              Hopping_Matrix(nf)%List(nc,1) = 1
              Hopping_Matrix(nf)%List(nc,2) = 1
              Hopping_Matrix(nf)%List(nc,3) = 1
              Hopping_Matrix(nf)%List(nc,4) = 0
           endif
           
           If (abs(Ham_Tperp) > Zero ) Then
              nc = nc + 1
              Hopping_Matrix(nf)%T(nc)    = cmplx(-Ham_Tperp,0.d0,kind(0.d0))
              Hopping_Matrix(nf)%List(nc,1) = 1
              Hopping_Matrix(nf)%List(nc,2) = 2
              Hopping_Matrix(nf)%List(nc,3) = 0
              Hopping_Matrix(nf)%List(nc,4) = 0
           endif

           If (abs(Ham_T2) > Zero ) Then
              nc = nc + 1
              Hopping_Matrix(nf)%T(nc)    = cmplx(-Ham_T2,0.d0,kind(0.d0))
              Hopping_Matrix(nf)%List(nc,1) = 2
              Hopping_Matrix(nf)%List(nc,2) = 2
              Hopping_Matrix(nf)%List(nc,3) = 0
              Hopping_Matrix(nf)%List(nc,4) = 1
              
              nc = nc + 1
              Hopping_Matrix(nf)%T(nc)    = cmplx(-Ham_T2,0.d0,kind(0.d0))
              Hopping_Matrix(nf)%List(nc,1) = 2
              Hopping_Matrix(nf)%List(nc,2) = 2
              Hopping_Matrix(nf)%List(nc,3) = 1
              Hopping_Matrix(nf)%List(nc,4) = 0
           endif


           Allocate ( Hopping_Matrix(nf)%T_Loc(Latt_Unit%Norb) )
           do nc = 1,Latt_Unit%Norb
              Hopping_Matrix(nf)%T_Loc(nc)  = cmplx(-Ham_Chem,0.d0,kind(0.d0))
           enddo
           Hopping_Matrix(nf)%N_Phi =  N_Phi
           Hopping_Matrix(nf)%Phi_X =  Phi_X
           Hopping_Matrix(nf)%Phi_Y =  Phi_Y
           Hopping_Matrix(nf)%Bulk =   Bulk
        enddo

        If  (Checkerboard)  then
           N_Fam  = 4
           if (abs(Ham_Tperp) > Zero )  N_Fam=5
           
           Allocate (L_FAM(N_FAM),  Prop_Fam(N_FAM), Multiplicity(Latt_unit%Norb) )
           
           No_Shift = 0
           If (abs(Ham_Tperp) > Zero ) No_Shift=1
           
           If     ( abs(Ham_T2)   <  Zero  .and. abs(Ham_Tperp) < Zero)    then
              L_FAM  = Latt%N/2
              Allocate (List_Fam(N_Fam,Latt%N/2,2))
              Multiplicity = 4
           elseif ( abs(Ham_T2)   <  Zero  .and. abs(Ham_Tperp) > Zero)    then
              L_FAM    = Latt%N/2
              L_FAM(5) = Latt%N
              Allocate (List_Fam(N_Fam,Latt%N,2))
              Multiplicity(1) = 5
              Multiplicity(2) = 1
           elseif ( abs(Ham_T2)   >  Zero  .and. abs(Ham_Tperp) < Zero)    then
              L_FAM    = Latt%N
              Allocate (List_Fam(N_Fam,Latt%N,2))
              Multiplicity = 4
           elseif ( abs(Ham_T2)   >  Zero  .and. abs(Ham_Tperp) > Zero)    then
              L_FAM    = Latt%N
              Allocate (List_Fam(N_Fam,Latt%N,2))
              Multiplicity = 5
              No_Shift     = 1
           endif
           L_FAM  = 0
           do I = 1,Latt%N
              if ( mod(Latt%List(I,1) + Latt%List(I,2),2) == 0 ) then
                 Nf = 1
                 L_FAM(Nf) = L_FAM(Nf) + 1
                 List_Fam(Nf,L_FAM(Nf),1) = I ! Unit cell
                 List_Fam(Nf,L_FAM(Nf),2) = 1 ! The bond (See above)
                 If (Abs(Ham_T2) > Zero) then
                    L_FAM(Nf) = L_FAM(Nf) + 1
                    List_Fam(Nf,L_FAM(Nf),1) = I ! Unit cell
                    List_Fam(Nf,L_FAM(Nf),2) = 3 + No_Shift ! The bond (See above)
                 endif
                 Nf = 2
                 L_FAM(Nf) = L_FAM(Nf) + 1
                 List_Fam(Nf,L_FAM(Nf),1) = I
                 List_Fam(Nf,L_FAM(Nf),2) = 2 
                 If (Abs(Ham_T2) > Zero) then
                    L_FAM(Nf) = L_FAM(Nf) + 1
                    List_Fam(Nf,L_FAM(Nf),1) = I ! Unit cell
                    List_Fam(Nf,L_FAM(Nf),2) = 4 + No_Shift ! The bond (See above)
                 endif
              else
                 Nf = 3
                 L_FAM(Nf) = L_FAM(Nf) + 1
                 List_Fam(Nf,L_FAM(Nf),1) = I
                 List_Fam(Nf,L_FAM(Nf),2) = 1  
                 If (Abs(Ham_T2) > Zero) then
                    L_FAM(Nf) = L_FAM(Nf) + 1
                    List_Fam(Nf,L_FAM(Nf),1) = I ! Unit cell
                    List_Fam(Nf,L_FAM(Nf),2) = 3 + No_Shift ! The bond (See above)
                 endif
                 Nf = 4
                 L_FAM(Nf) = L_FAM(Nf) + 1
                 List_Fam(Nf,L_FAM(Nf),1) = I
                 List_Fam(Nf,L_FAM(Nf),2) = 2  
                 If (Abs(Ham_T2) > Zero) then
                    L_FAM(Nf) = L_FAM(Nf) + 1
                    List_Fam(Nf,L_FAM(Nf),1) = I ! Unit cell
                    List_Fam(Nf,L_FAM(Nf),2) = 4 + No_Shift ! The bond (See above)
                 endif
              endif
              If (Abs(Ham_Tperp) > Zero) then
                 Nf = 5
                 L_FAM(Nf) = L_FAM(Nf) + 1
                 List_Fam(Nf,L_FAM(Nf),1) = I
                 List_Fam(Nf,L_FAM(Nf),2) = 3  
              Endif
           enddo
           ! Test
           If (Test) then
              Write(6,*)  N_FAM,  L_FAM
              Write(6,*)  Ham_T1,Ham_T2, Ham_Tperp
              Do nf = 1,N_FAM
                 Do n = 1,L_Fam(nf)
                    I =  List_Fam(Nf,n,1)
                    nb = List_Fam(Nf,n,2)
                    Write(6,"(I3,2x,I3,2x,I3,2x,I3,2x,I3,2x,I3)")   Latt%list(I,1), Latt%list(I,2), Hopping_Matrix(1)%List(nb,1),Hopping_Matrix(1)%List(nb,2), &
                         &Hopping_Matrix(1)%List(nb,3), Hopping_Matrix(1)%List(nb,4)
                 enddo
                 Write(6,*)
              enddo
           endif
        endif


        
      end Subroutine Set_Default_hopping_parameters_Bilayer_square
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>
!> @brief 
!> Given the checkerbord decompostion  
!> the routine allocates and sets OP_T Set_Default_hopping_parameters_"Lattice" routine,
!> this routine generates the data for the symmetric decompostion. 
!     
!--------------------------------------------------------------------      
      Subroutine Symmetrize_Families
        implicit none
        
        ! In Families.  Out Symmetrized Families.

        !  Make a copy  of the unsymmetrized forms
        Integer                              ::  N_Fam_C
        Integer, allocatable                 ::  L_Fam_C(:),  List_Fam_C(:,:,:)
        Real (Kind=Kind(0.d0)), allocatable  ::  Prop_Fam_C(:)

        Integer :: n,n1,n2, n_f_max, n_l_max, nc
        Integer, allocatable ::  list_Fam_tmp(:)
        
        ! Copy
        N_FAM_C = N_FAM
        Allocate(L_FAM_C(N_FAM_C)) 
        n2 = Size(List_Fam,2)
        Allocate ( List_Fam_C(N_FAM_C,n2,2), Prop_Fam_C(N_Fam_C) )
        L_FAM_C    = L_FAM
        List_Fam_C = List_Fam
        Prop_Fam_C = Prop_Fam

        ! Re-allocate
        N_FAM  =  2*N_FAM_C - 1
        Deallocate (L_Fam, List_Fam, Prop_Fam)
        Allocate (L_Fam(N_Fam), List_Fam(N_FAM,n2,2), Prop_Fam(N_Fam) )

        ! Symmetrize
        ! Find de longest family.
        n_l_max = 0
        do n = 1, N_FAM_C
           if (L_FAM_C(n) > n_l_max ) then
              n_l_max = L_FAM_C(n)
              n_f_max = n
           endif
        enddo

        Allocate( list_Fam_tmp(N_FAM) )
        nc = 0
        Do n = 1, N_FAM_C   
           nc = nc + 1
           list_Fam_tmp(nc) = n
        Enddo
        Do n = N_FAM_C-1,1,-1 
           nc = nc + 1
           list_Fam_tmp(nc) = n
        Enddo
        ! Place the largest familly in the middle and set the time step.
        Prop_Fam          = 0.5D0
        Prop_Fam(N_Fam_C) = 1.D0
        list_Fam_tmp(N_FAM_C) = n_f_max
        list_Fam_tmp(1)       = N_Fam_C
        list_Fam_tmp(N_FAM )  = N_Fam_C

        do n = 1,N_FAM
           n1 = list_Fam_tmp(n)
           L_FAM(n)        = L_FAM_C(n1)
           List_Fam(n,:,:) = List_Fam_C(n1,:,:)
        enddo

        ! Clean
        Deallocate( L_FAM_C, List_Fam_C, Prop_Fam_C, List_Fam_tmp )

        Write(6,*)  N_FAM
        Write(6,*)  L_FAM
        Write(6,*)  Prop_Fam
        
      end Subroutine Symmetrize_Families
      
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>
!> @brief 
!> Given the Hopping-martix, and if required the checkerboard decomposion  (i.e. private data of this module) 
!> the routine allocates and sets OP_T
!
!--------------------------------------------------------------------      
      Subroutine Predefined_Hopping_new(List,Invlist,Latt,  Latt_unit,  Dtau,Checkerboard, Symm,  OP_T )
        
        Implicit none
        
        Integer, Intent(IN), Dimension(:,:)                 :: List, Invlist
        Type(Lattice),  Intent(in)                          :: Latt
        Type(Unit_cell),Intent(in)                          :: Latt_unit
        Real (Kind=Kind(0.d0)), Intent(In)                  :: Dtau
        Logical                                             :: Checkerboard, Symm
        
        Type(Operator), Intent(Out),  dimension(:,:), allocatable  :: Op_T 
        


        ! Local 
        Integer                           :: Ndim, N_FL, N_Phi, I, J, I1, J1, no_I, no_J, nf
        Integer                           :: n_1, n_2, Nb, n_f,l_f, n_l, N, nc
        Real   (Kind=Kind(0.d0))          :: Ham_T, Ham_Chem,  Phi_X, Phi_Y
        Logical                           :: Bulk
        Complex(Kind=Kind(0.d0))          :: Z 
        
        N_FL =  size(Hopping_Matrix,1)
        !Write(6,*)  'N_FL ', N_FL
        Ndim =  Latt%N * Latt_Unit%Norb
        
        If ( .not. Checkerboard) then
           allocate(Op_T(1,N_FL))
           do nf = 1,N_FL
              Call Op_make(Op_T(1,nf),Ndim)   ! This is too restrictive for the  Kondo type models. The hopping only occurs on one  subsystem. 
              N_Phi     = Hopping_Matrix(nf)%N_Phi
              Phi_X     = Hopping_Matrix(nf)%Phi_X  
              Phi_Y     = Hopping_Matrix(nf)%Phi_Y 
              Bulk      = Hopping_Matrix(nf)%Bulk 
              DO I = 1, Latt%N
                 do Nb = 1, Hopping_Matrix(nf)%N_bonds
                    no_I = Hopping_Matrix(nf)%list(Nb,1)
                    no_J = Hopping_Matrix(nf)%list(Nb,2)
                    n_1  = Hopping_Matrix(nf)%list(Nb,3)
                    n_2  = Hopping_Matrix(nf)%list(Nb,4)
                    J    = Latt%nnlist(I,n_1,n_2)
                    Z    = Generic_hopping(I,no_I, n_1, n_2, no_J, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit)
                    I1   = Invlist(I,no_I)
                    J1   = Invlist(J,no_J)
                    Op_T(1,nf)%O(I1,J1) = Hopping_Matrix(nf)%T(Nb)*Z
                    Op_T(1,nf)%O(J1,I1) = Conjg(Hopping_Matrix(nf)%T(Nb)*Z)
                 enddo
                 ! T(N_b=1..N_bonds)
                 ! List(N_b,1) = no_1
                 ! List(N_b,2) = no_2
                 ! List(N_b,3) = n_1
                 ! List(N_b,4) = n_2
                 ! H_[(i,no_1),(i + n_1 a_1 + n_2 a_2,no_2)] = T(N_b)
                 Do no_I = 1, Latt_Unit%Norb
                    I1   = Invlist(I,no_I)
                    Op_T(1,nf)%O(I1,I1) = Hopping_Matrix(nf)%T_Loc(no_I)
                 Enddo
              enddo
              Do I = 1,Ndim
                 Op_T(1,nf)%P(i) = i 
              Enddo
              if ( Hopping_Matrix(nf)%N_bonds == 0 ) then 
                 Op_T(1,nf)%g = 0.d0
              else
                 Op_T(1,nf)%g = -Dtau
              endif
              Op_T(1,nf)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
              Call Op_set(Op_T(1,nf))
              !Do I = 1,Size(Op_T(nc,n)%E,1)
              !   Write(6,*) Op_T(nc,n)%E(I)
              !Enddo
           Enddo
        Elseif (Checkerboard) then
           If (Symm) Call Symmetrize_families
           N = 0
           do n_f = 1,N_Fam
              N = N +  L_Fam(n_f)
           enddo
           allocate(Op_T(N,N_FL))
           do nf = 1,N_FL
              N_Phi     = Hopping_Matrix(nf)%N_Phi
              Phi_X     = Hopping_Matrix(nf)%Phi_X  
              Phi_Y     = Hopping_Matrix(nf)%Phi_Y 
              Bulk      = Hopping_Matrix(nf)%Bulk 
              do nc = 1, Size(Op_T,1)
                 Call Op_make(Op_T(nc,nf),2)
              enddo
              nc = 0
              Do n_f = 1, N_Fam
                 Do l_f = 1, L_Fam(n_f)
                    I  = List_Fam(n_f,l_f,1)
                    nb = List_Fam(n_f,l_f,2)
                    no_I = Hopping_Matrix(nf)%list(Nb,1)
                    no_J = Hopping_Matrix(nf)%list(Nb,2)
                    n_1  = Hopping_Matrix(nf)%list(Nb,3)
                    n_2  = Hopping_Matrix(nf)%list(Nb,4)
                    J    = Latt%nnlist(I,n_1,n_2)
                    Z    = Generic_hopping(I,no_I, n_1, n_2, no_J, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit)
                    I1   = Invlist(I,no_I)
                    J1   = Invlist(J,no_J)
                    nc = nc + 1
                    Op_T(nc,nf)%P(1) = I1 
                    Op_T(nc,nf)%P(2) = J1 
                    Op_T(nc,nf)%O(1,2) = Hopping_Matrix(nf)%T(Nb)*Z
                    Op_T(nc,nf)%O(2,1) = Conjg(Hopping_Matrix(nf)%T(Nb)*Z)
                    Op_T(nc,nf)%O(1,1) = Hopping_Matrix(nf)%T_loc(no_I)/Multiplicity(no_I)
                    Op_T(nc,nf)%O(2,2) = Hopping_Matrix(nf)%T_loc(no_J)/Multiplicity(no_J)
                    Op_T(nc,nf)%g = -Dtau*Prop_Fam(n_f)
                    Op_T(nc,nf)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
                    Call Op_set(Op_T(nc,nf))
                 enddo
              enddo
           enddo
        else
           Write(6,*) 'Generic checkerboard  is not yet implemented'
           Stop
        endif
        
      end Subroutine Predefined_Hopping_new

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>
!> @brief 
!> The subroutine computes the kinietic energy  based on the generic form of the 
!> the hopping matrix.
!> 
!--------------------------------------------------------------------
      Subroutine  Predefined_Hop_Compute_Kin(List,Invlist, Latt, Latt_unit, GRC, Z_Kin)

        Implicit none

        Integer, Intent(IN), Dimension(:,:)                 :: List, Invlist
        Type(Lattice),  Intent(in)                          :: Latt
        Type(Unit_cell),Intent(in)                          :: Latt_unit
        Complex (Kind=Kind(0.d0)), intent(in), Dimension(:,:,:) :: GRC(:,:,:)
        Complex (Kind=Kind(0.d0)),  intent(out) :: Z_kin

        !Local
        Integer                           :: Ndim, N_FL, N_Phi, I, J, I1, J1, no_I, no_J, nf
        Integer                           :: n_1, n_2, Nb, n_f,l_f, n_l, N, nc
        Real   (Kind=Kind(0.d0))          :: Ham_T, Ham_Chem,  Phi_X, Phi_Y
        Logical                           :: Bulk
        Complex(Kind=Kind(0.d0))          :: Z 

        N_FL  =  Size(GRC,3)

        Z_Kin = cmplx(0.d0,0.d0,Kind(0.d0))
        do nf = 1,N_FL
           N_Phi     = Hopping_Matrix(nf)%N_Phi
           Phi_X     = Hopping_Matrix(nf)%Phi_X  
           Phi_Y     = Hopping_Matrix(nf)%Phi_Y 
           Bulk      = Hopping_Matrix(nf)%Bulk 
           DO I = 1, Latt%N
              do Nb = 1, Hopping_Matrix(nf)%N_bonds
                 no_I = Hopping_Matrix(nf)%list(Nb,1)
                 no_J = Hopping_Matrix(nf)%list(Nb,2)
                 n_1  = Hopping_Matrix(nf)%list(Nb,3)
                 n_2  = Hopping_Matrix(nf)%list(Nb,4)
                 J    = Latt%nnlist(I,n_1,n_2)
                 Z    = Generic_hopping(I,no_I, n_1, n_2, no_J, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit)
                 I1   = Invlist(I,no_I)
                 J1   = Invlist(J,no_J)
                 Z_Kin = Z_Kin + Hopping_Matrix(nf)%T(Nb)*Z * GRC(I1,J1,nf) + conjg(Hopping_Matrix(nf)%T(Nb)*Z)*GRC(J1,I1,nf)
              enddo
              Do no_I = 1, Latt_Unit%Norb
                 I1   = Invlist(I,no_I)
                 Z_Kin = Z_Kin   +  Hopping_Matrix(nf)%T_Loc(no_I)*GRC(I1,I1,nf)
              Enddo
           enddo
        enddo
           
      end Subroutine Predefined_Hop_Compute_Kin
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>
!> @brief 
!>    Hopping, with ot without checkerboard  
!>    Per flavor, the  hopping is given by
!>    \f[  e^{ - \Delta \tau  H_t  }   = \prod_{n=1}^{N} e^{ - \Delta \tau_n  H_t(n) }   \f]
!>    If  _Symm_ is set to true and if  _Checkeborad_ is on, then  one will carry out a
!>    symmetric decomposition so as to preserve  the hermitian properties of the hopping.
!>    Thereby   OP_T has dimension OP_T(N,N_FL)
!> @param [in]  Latttice_type
!>    Character(64)
!>\verbatim 
!>     Square,  Honeycomb, Pi_Flux 
!>\endverbatim 
!> @param [in]  Latt_unit
!>    Type(Unit_cell)
!> \verbatim
!>     Contains number of orbitals per unit cell and positions, as well as coordination number
!> \endverbatim
!> @param [in]  Ndim
!>    Integer
!> \verbatim
!>     Number of orbitals      
!> \endverbatim
!> @param [in]  List, Invlist
!>    Integer(:,:)
!> \verbatim
!>      List(I=1.. Ndim,1)    =   Unit cell of site I    
!>      List(I=1.. Ndim,2)    =   Orbital index  of site I    
!>      Invlist(Unit_cell,Orbital) = site I    
!> \endverbatim
!> @param [in]    Latt
!>    Type(Lattice)
!> \verbatim
!>      The Lattice
!> \endverbatim
!> @param [in]  Dtau
!>    Real
!> \verbatim
!>      Imaginary time step
!> \endverbatim
!> @param [in]  Ham_T
!>    Real
!> \verbatim
!>      Hopping matrix element
!> \endverbatim
!> @param [in]  Ham_Chem
!>    Real
!> \verbatim
!>      Chemical potential
!> \endverbatim
!> @param [in]  XB_X, YB_Y
!>    Real
!> \verbatim
!>      X, Y  Boundary conditions
!> \endverbatim
!> @param [in]  Phi_X, Phi_Y 
!>    Real
!> \verbatim
!>      X, Y  Fluxes
!> \endverbatim
!> @param [in]  N_FL
!>    Integer
!> \verbatim
!>      Flavors
!> \endverbatim
!> @param [in]  Checkerboard
!>    Logical
!> \verbatim
!>      Allows for checkerboard decomposition
!> \endverbatim
!> @param [in]  Symm
!>    Logical
!> \verbatim
!>      Allows for symmetric checkerboard decomposition
!> \endverbatim
!> @param [out]  OP_T 
!>    Type(operator)(N,N_FL)
!> \verbatim
!>      Hopping
!> \endverbatim
!> @param [in]  Dimer 
!>    Real, Optional.  Modulation of hopping that breaks lattice symmetries so as to generate a unique
!>    ground state for the half-filled case.  This option is  effective only  if the checkerboard
!>    decomposition is not used. It is presently implemented for the square and one-dimensional lattices.
!> \verbatim
!>      Hopping
!> \endverbatim
!>       
!------------------------------------------------------------------
      Subroutine Predefined_Hopping(Lattice_type, Ndim, List,Invlist,Latt,  Latt_unit,  &
           &                        Dtau, Ham_T, Ham_Chem, Phi_X, Phi_Y, Bulk,  N_Phi, &
           &                        N_FL,  Checkerboard, Symm, OP_T, Dimer )

        Implicit none

        Character (len=64), Intent(IN)               :: Lattice_type
        Integer, Intent(IN)                          :: Ndim, N_FL, N_Phi
        Integer, Intent(IN), Dimension(:,:)          :: List, Invlist
        Type(Lattice),  Intent(in)                   :: Latt
        Type(Unit_cell),Intent(in)                   :: Latt_unit
        Real (Kind=Kind(0.d0)), Intent(In)           :: Dtau, Ham_T, Ham_Chem,  Phi_X, Phi_Y
        Logical                                      :: Checkerboard, Symm,  Bulk
        Real(Kind=Kind(0.d0)), Intent(IN), Optional  :: Dimer
        
        Type(Operator), Intent(Out),  dimension(:,:), allocatable  :: Op_T 


        !Local
        Integer :: I, I1, J1, I2, n, Ncheck,nc, nc1, no, N_Fam, L_FAM, Ix, Iy, n1
        Complex (Kind=Kind(0.d0)) :: ZX, ZY, Z
        Real    (Kind=Kind(0.d0)) :: del_p(2), X, g
        
        If ( .not. Checkerboard) then
           allocate(Op_T(1,N_FL))
           do n = 1,N_FL
              Call Op_make(Op_T(1,n),Ndim)
              nc = 1
              Select case (Lattice_type)
              Case ("Square")
                 If ( Latt_unit%N_coord == 2 ) then   !  This is for the 2D case
                    DO I = 1, Latt%N
                       I1 = Latt%nnlist(I,1,0)
                       I2 = Latt%nnlist(I,0,1)
                       ZX = Generic_hopping(I,1, 1, 0, 1, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                       ZY = Generic_hopping(i,1, 0, 1, 1, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit)
                       Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))*ZX
                       Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T, 0.d0, kind(0.D0))*conjg(ZX)
                       Op_T(nc,n)%O(I,I2) = cmplx(-Ham_T, 0.d0, kind(0.D0))*ZY
                       Op_T(nc,n)%O(I2,I) = cmplx(-Ham_T, 0.d0, kind(0.D0))*conjg(ZY)
                       Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                    Enddo
                    if (Present(Dimer)) then
                       Do I = 1, Latt%N
                          Ix = Latt%list(I,1); Iy = Latt%list(I,2)
                          If ( mod(Ix + Iy,2)  == 0 ) then
                             I1 = I
                             I2 = Latt%nnlist(I,1,0)
                             Op_T(nc,n)%O(I1,I2)  = Op_T(nc,n)%O(I1,I2) + Dimer
                             Op_T(nc,n)%O(I2,I1)  = Op_T(nc,n)%O(I2,I1) + Dimer
                          endif
                       enddo
                    endif
                 else  ! One dimensional
                    DO I = 1, Latt%N
                       I1 = Latt%nnlist(I,1,0)
                       ZX = Generic_hopping(I,1, 1, 0, 1, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                       Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))*ZX
                       Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T, 0.d0, kind(0.D0))*conjg(ZX)
                       Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                    Enddo
                    if (Present(Dimer)) then
                       Do I = 1, Latt%N
                          Ix = Latt%list(I,1); Iy = Latt%list(I,2)
                          If ( mod(Ix + Iy,2)  == 0 ) then
                             I1 = I
                             I2 = Latt%nnlist(I,1,0)
                             Op_T(nc,n)%O(I1,I2)  = Op_T(nc,n)%O(I1,I2) + Dimer
                             Op_T(nc,n)%O(I2,I1)  = Op_T(nc,n)%O(I2,I1) + Dimer
                          endif
                       enddo
                    endif
                 endif
              Case ("Honeycomb")
                 DO I = 1, Latt%N
                    do no = 1,Latt_unit%Norb
                       I1 = Invlist(I,no)
                       Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                    enddo
                    I1 = Invlist(I,1)
                    J1 = I1
                    Do nc1 = 1,Latt_unit%N_coord
                       select case (nc1)
                       case (1)
                          J1 = invlist(I,2)
                          ZX = Generic_hopping(I,1, 0, 0, 2, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                       case (2)
                          J1 = invlist(Latt%nnlist(I,1,-1),2)
                          ZX = Generic_hopping(I,1, 1, -1, 2, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                       case (3)
                          J1 = invlist(Latt%nnlist(I,0,-1),2) 
                          ZX = Generic_hopping(I,1, 0, -1, 2, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                       case default
                          Write(6,*) ' Error in  Ham_Hop '  
                          Stop
                       end select
                       Op_T(nc,n)%O(I1,J1) = cmplx(-Ham_T,    0.d0, kind(0.D0)) * ZX
                       Op_T(nc,n)%O(J1,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0)) * CONJG(ZX) 
                    Enddo
                 Enddo
              case("Pi_Flux")
                 DO I = 1, Latt%N
                    do no = 1,Latt_unit%Norb
                       I1 = Invlist(I,no)
                       Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                    enddo
                    I1 = Invlist(I,1)
                    J1 = I1
                    Do nc1 = 1,Latt_unit%N_coord
                       select case (nc1)
                       case (1)
                          J1 = invlist(I,2) 
                       case (2)
                          J1 = invlist(Latt%nnlist(I,0, 1),2) 
                       case (3)
                          J1 = invlist(Latt%nnlist(I,-1,1),2) 
                       case (4)
                          J1 = invlist(Latt%nnlist(I,-1,0),2) 
                       case default
                          Write(6,*) ' Error in  Ham_Hop '  
                          Stop
                       end select
                       if (nc1 == 1 ) then
                          Op_T(nc,n)%O(I1,J1) = cmplx( Ham_T,    0.d0, kind(0.D0))
                          Op_T(nc,n)%O(J1,I1) = cmplx( Ham_T,    0.d0, kind(0.D0))
                       Else
                          Op_T(nc,n)%O(I1,J1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                          Op_T(nc,n)%O(J1,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                       endif
                    Enddo
                 Enddo
              case default 
                 Write(6,*) "Lattice not yet implemented!"
                 Stop
              end select
              Do I = 1,Ndim
                 Op_T(nc,n)%P(i) = i 
              Enddo
              if ( abs(Ham_T) < 1.E-6  .and.  abs(Ham_chem) < 1.E-6 ) then 
                 Op_T(nc,n)%g = 0.d0
              else
                 Op_T(nc,n)%g = -Dtau
              endif
              Op_T(nc,n)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
              Call Op_set(Op_T(nc,n))
              !Do I = 1,Size(Op_T(nc,n)%E,1)
              !   Write(6,*) Op_T(nc,n)%E(I)
              !Enddo
           Enddo
        else
           select case (Lattice_type)
           case ("Square")
              If (Symm) then
                 N_Fam = 2*(2*Latt_Unit%N_coord) -1  !  Number of Families = 7
                 L_Fam = Latt%N/2                    !  Length of Families = LQ/2
                 Allocate(Op_T(N_FAM*L_FAM,N_FL))
                 do n = 1,N_FL
                    do I  =  1, N_FAM*L_FAM
                       call Op_make(Op_T(I,n),2)
                    enddo
                    nc = 0
                    do nc1 = 1,N_Fam
                       do I = 1,Latt%N
                          if ( mod(Latt%List(I,1) + Latt%List(I,2),2) == 0 ) then
                             I1 = I
                             nc = nc + 1
                             select case (nc1)
                             case(1)
                                I2 = latt%nnlist(I1, 1, 0) 
                                ZX = Generic_hopping(I,1,  1,  0, 1, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                                g = -Dtau/2.d0
                             case(2)
                                I2 = latt%nnlist(I1, 0, 1)
                                ZX = Generic_hopping(I,1,  0,  1, 1, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                                g = -Dtau/2.d0
                             case(3)
                                I2 = latt%nnlist(I1,-1, 0)
                                ZX = Generic_hopping(I,1, -1,  0, 1, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                                g = -Dtau/2.d0
                             case(4)
                                I2 = latt%nnlist(I1, 0,-1)
                                ZX = Generic_hopping(I,1,  0, -1, 1, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                                g = -Dtau
                             case(5)
                                I2 = latt%nnlist(I1,-1, 0)
                                ZX = Generic_hopping(I,1, -1,  0, 1, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                                g = -Dtau/2.d0
                             case(6)
                                I2 = latt%nnlist(I1, 0, 1)
                                ZX = Generic_hopping(I,1,  0,  1, 1, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                                g = -Dtau/2.d0
                             case(7)
                                I2 = latt%nnlist(I1, 1, 0) 
                                ZX = Generic_hopping(I,1,  1,  0, 1, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                                g = -Dtau/2.d0
                             end select
                             !Write(6,*) nc,nc1, Latt%List(I1,1), Latt%List(I1,2),Latt%List(I2,1), Latt%List(I2,2), I1, I2
                             Op_T(nc,n)%P(1)   = I1
                             Op_T(nc,n)%P(2)   = I2
                             Op_T(nc,n)%O(1,2) = cmplx(-Ham_T ,0.d0, kind(0.D0))*ZX 
                             Op_T(nc,n)%O(2,1) = cmplx(-Ham_T ,0.d0, kind(0.D0))*CONJG(ZX)
                             Op_T(nc,n)%O(1,1) = cmplx(-Ham_Chem/4.d0 ,0.d0, kind(0.D0)) 
                             Op_T(nc,n)%O(2,2) = cmplx(-Ham_Chem/4.d0 ,0.d0, kind(0.D0))
                             if ( abs(Ham_T) < 1.E-6  .and.  abs(Ham_chem) < 1.E-6 ) then 
                                Op_T(nc,n)%g = 0.d0
                             else
                                Op_T(nc,n)%g = g
                             endif
                             Op_T(nc,n)%alpha  = cmplx( 0.d0, 0.d0, kind(0.D0) )
                             Call Op_set( Op_T(nc,n) )
                          endif
                       Enddo
                    Enddo
                 Enddo
              Else   !  Trotter no symm
                 N_Fam = 2*Latt_unit%N_coord  !  Number of Families = 4
                 L_Fam = Latt%N/2             !  Length of Families = LQ/2
                 Allocate(Op_T(N_FAM*L_FAM,N_FL))
                 do n = 1,N_FL
                    do I  =  1, N_FAM*L_FAM
                       call Op_make(Op_T(I,n),2)
                    enddo
                    nc = 0
                    do nc1 = 1,N_Fam
                       do I = 1,Latt%N
                          if ( mod(Latt%List(I,1) + Latt%List(I,2),2) == 0 ) then
                             I1 = I
                             nc = nc + 1
                             select case (nc1)
                             case(1)
                                I2 = latt%nnlist(I1, 1, 0) 
                                ZX = Generic_hopping(I,1,  1,  0, 1, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                             case(2)
                                I2 = latt%nnlist(I1, 0, 1)
                                ZX = Generic_hopping(I,1,  0,  1, 1, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                             case(3)
                                I2 = latt%nnlist(I1,-1, 0)
                                ZX = Generic_hopping(I,1, -1,  0, 1, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                             case(4)
                                I2 = latt%nnlist(I1, 0,-1)
                                ZX = Generic_hopping(I,1,  0, -1, 1, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                             end select
                             !Write(6,*) nc,nc1, Latt%List(I1,1), Latt%List(I1,2),Latt%List(I2,1), Latt%List(I2,2), I1, I2
                             Op_T(nc,n)%P(1) = I1
                             Op_T(nc,n)%P(2) = I2
                             Op_T(nc,n)%O(1,2) = cmplx(-Ham_T ,0.d0, kind(0.D0)) * ZX
                             Op_T(nc,n)%O(2,1) = cmplx(-Ham_T ,0.d0, kind(0.D0)) * CONJG(ZX)
                             Op_T(nc,n)%O(1,1) = cmplx(-Ham_Chem/4.d0 ,0.d0, kind(0.D0)) 
                             Op_T(nc,n)%O(2,2) = cmplx(-Ham_Chem/4.d0 ,0.d0, kind(0.D0))
                             if ( abs(Ham_T) < 1.E-6  .and.  abs(Ham_chem) < 1.E-6 ) then 
                                Op_T(nc,n)%g = 0.d0
                             else
                                Op_T(nc,n)%g = -Dtau
                             endif
                             Op_T(nc,n)%alpha  = cmplx( 0.d0, 0.d0, kind(0.D0) )
                             Call Op_set( Op_T(nc,n) )
                          endif
                       Enddo
                    Enddo
                 Enddo
              Endif
           case ("Pi_Flux")
              If (Symm) Then
                 Write(6,*) 'Pi-Flux Symm is not yet implemented'
                 Stop
              else
                 Allocate(Op_T(Latt_unit%N_coord*Latt%N,N_FL))
                 do n = 1,N_FL
                    !Write(6,*) 'N_coord, Latt%N ',  N_coord, Latt%N
                    do i  =  1, Latt_unit%N_coord*Latt%N
                       call Op_make(Op_T(i,n),2)
                    enddo
                    nc = 0
                    do nc1 = 1,Latt_unit%N_coord
                       do I = 1,Latt%N
                          I1 = Invlist(I,1)
                          nc = nc + 1
                          If (nc1 == 1 )  I2 = invlist(I,2)
                          If (nc1 == 2 )  I2 = invlist(Latt%nnlist(I,0, 1),2) 
                          If (nc1 == 3 )  I2 = invlist(Latt%nnlist(I,-1,1),2)
                          If (nc1 == 4 )  I2 = invlist(Latt%nnlist(I,-1,0),2) 
                          Op_T(nc,n)%P(1) = I1
                          Op_T(nc,n)%P(2) = I2
                          if (nc1 == 1 ) then
                             Op_T(nc,n)%O(1,2) = cmplx( Ham_T,    0.d0, kind(0.D0))
                             Op_T(nc,n)%O(2,1) = cmplx( Ham_T,    0.d0, kind(0.D0))
                          Else
                             Op_T(nc,n)%O(1,2) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                             Op_T(nc,n)%O(2,1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                          endif
                          Op_T(nc,n)%O(1,1) = cmplx(-Ham_Chem/4.d0 ,0.d0, kind(0.D0)) 
                          Op_T(nc,n)%O(2,2) = cmplx(-Ham_Chem/4.d0 ,0.d0, kind(0.D0))
                          if ( abs(Ham_T) < 1.E-6  .and.  abs(Ham_chem) < 1.E-6 ) then 
                             Op_T(nc,n)%g = 0.d0
                          else
                             Op_T(nc,n)%g = -Dtau
                          endif
                          Op_T(nc,n)%alpha  = cmplx( 0.d0, 0.d0, kind(0.D0) )
                          Call Op_set( Op_T(nc,n) )
                       Enddo
                    Enddo
                 Enddo
              Endif
           case ("Honeycomb")
              If (Symm) then
                 Allocate(Op_T((2*Latt_unit%N_coord-1)*Latt%N,N_FL))
                 do n = 1,N_FL
                    do i  =  1, (2*Latt_unit%N_coord-1)*Latt%N
                       call Op_make(Op_T(i,n),2)
                    enddo
                    nc = 0
                    do nc1 = 1,(2*Latt_unit%N_coord -1) 
                       do I = 1,Latt%N
                          I1 = invlist(I,1)
                          nc = nc + 1
                          select case (nc1)
                          case(1)
                             I2 = invlist(I,2)
                             ZX = Generic_hopping(I,1, 0, 0, 2, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                             g = -Dtau/2.d0
                          case(2)
                             I2 = invlist(latt%nnlist(I,1,-1),2)
                             ZX = Generic_hopping(I,1, 1, -1, 2, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                             g = -Dtau/2.d0
                          case(3)
                             I2 = invlist(latt%nnlist(I,0,-1),2)
                             ZX = Generic_hopping(I,1, 0, -1, 2, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                             g = -Dtau
                          case(4)
                             I2 = invlist(latt%nnlist(I,1,-1),2)
                             ZX = Generic_hopping(I,1, 1, -1, 2, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                             g = -Dtau/2.d0
                          case(5)
                             I2 = invlist(I,2)
                             ZX = Generic_hopping(I,1, 0,  0, 2, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                             g = -Dtau/2.d0
                          end select
                          Op_T(nc,n)%P(1) = I1
                          Op_T(nc,n)%P(2) = I2
                          Op_T(nc,n)%O(1,2) = cmplx(-Ham_T ,0.d0, kind(0.D0)) * ZX
                          Op_T(nc,n)%O(2,1) = cmplx(-Ham_T ,0.d0, kind(0.D0)) * CONJG(ZX)
                          Op_T(nc,n)%O(1,1) = cmplx(-Ham_Chem/3.d0 ,0.d0, kind(0.D0)) 
                          Op_T(nc,n)%O(2,2) = cmplx(-Ham_Chem/3.d0 ,0.d0, kind(0.D0))
                          if ( abs(Ham_T) < 1.E-6  .and.  abs(Ham_chem) < 1.E-6 ) then 
                             Op_T(nc,n)%g = 0.d0
                          else
                             Op_T(nc,n)%g = g
                          endif
                          Op_T(nc,n)%alpha  = cmplx( 0.d0, 0.d0, kind(0.D0) )
                          Call Op_set( Op_T(nc,n) )
                       Enddo
                    Enddo
                 Enddo
              Else
                 Allocate(Op_T(Latt_unit%N_coord*Latt%N,N_FL))
                 do n = 1,N_FL
                    do i  =  1, Latt_unit%N_coord*Latt%N
                       call Op_make(Op_T(i,n),2)
                    enddo
                    nc = 0
                    do nc1 = 1,Latt_unit%N_coord
                       do I = 1,Latt%N
                          I1 = invlist(I,1)
                          nc = nc + 1
                          select case (nc1)
                          case(1)
                             I2 = invlist(I,2)
                             ZX = Generic_hopping(I,1, 0, 0, 2, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                          case(2)
                             I2 = invlist(latt%nnlist(I,1,-1),2)
                             ZX = Generic_hopping(I,1, 1, -1, 2, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                          case(3)
                             I2 = invlist(latt%nnlist(I,0,-1),2)
                             ZX = Generic_hopping(I,1, 0, -1, 2, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit) 
                          end select
                          Op_T(nc,n)%P(1) = I1
                          Op_T(nc,n)%P(2) = I2
                          Op_T(nc,n)%O(1,2) = cmplx(-Ham_T ,0.d0, kind(0.D0)) * ZX 
                          Op_T(nc,n)%O(2,1) = cmplx(-Ham_T ,0.d0, kind(0.D0)) * CONJG(ZX)
                          Op_T(nc,n)%O(1,1) = cmplx(-Ham_Chem/3.d0 ,0.d0, kind(0.D0)) 
                          Op_T(nc,n)%O(2,2) = cmplx(-Ham_Chem/3.d0 ,0.d0, kind(0.D0))
                          if ( abs(Ham_T) < 1.E-6  .and.  abs(Ham_chem) < 1.E-6 ) then 
                             Op_T(nc,n)%g = 0.d0
                          else
                             Op_T(nc,n)%g = -Dtau
                          endif
                          Op_T(nc,n)%alpha  = cmplx( 0.d0, 0.d0, kind(0.D0) )
                          Call Op_set( Op_T(nc,n) )
                       Enddo
                    Enddo
                 Enddo
              Endif
              case default
                 Write(6,*) "Checkeboard is not implemented for this lattice"
              Stop
           End select
        endif
        
      End Subroutine Predefined_Hopping

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>
!> @brief 
!> This  function provides generic hopping.
!>    
!--------------------------------------------------------------------

      complex  (Kind=kind(0.d0)) function Generic_hopping(i,no_i, del_1, del_2, no_j, N_Phi, Flux_1,Flux_2, Bulk, Latt, Latt_Unit)

        Use Lattices_v3 
        Implicit none

        
        Integer        ,  Intent(In) :: N_Phi, i, no_i, del_1, del_2, no_j
        Type(Unit_cell),  Intent(In) :: Latt_Unit
        Type(Lattice)  ,  Intent(In) :: Latt
        Real (Kind = Kind(0.d0)), intent(In) :: Flux_1,Flux_2
        Logical        ,  Intent(In) :: Bulk 


        !Local
        Integer                   :: j, N1, N2
        real (Kind=Kind(0.d0))    :: xj_p(2), xi_p(2), xjp_p(2), del_p(2), A_p(2), pi, XB_p(2), V, B, Zero, x_p(2), x1_p(2)

        Complex (Kind=Kind(0.d0)) :: Z_hop

        
        Z_hop = cmplx(1.d0,0.d0,kind(0.d0))
        
        xj_p =  real(latt%list(i,1) + del_1 ,kind(0.d0)) * latt%a1_p  +  real(latt%list(i,2) + del_2 ,kind(0.d0)) * latt%a2_p
        ! Check if you have crossed the boundary:  xj_p  = xjp_p + N1*L1_p  + N2*L2_p  with  xjp_p  in the set of lattice sites. 
        N1 = 0; N2 = 0
        Call npbc(xjp_p, xj_p, Latt%L1_p, Latt%L2_p,  N1, N2)
        XB_p = real(N1,kind(0.d0))*Latt%L1_p  +  real(N2,kind(0.d0))*Latt%L2_p  
        xj_p (:) = xj_p (:) + Latt_unit%Orb_pos_p(no_j,:)
        xjp_p(:) = xjp_p(:) + Latt_unit%Orb_pos_p(no_j,:)   
        xi_p    = real(latt%list(i,1), kind(0.d0)) * latt%a1_p  +  real(latt%list(i,2),kind(0.d0)) * latt%a2_p
        xi_p(:) = xi_p(:) +  Latt_unit%Orb_pos_p(no_i,:)

        !!Check that  xjp_p(:) + XB_p(:) =  xj_p(:)
        !!x1_p(:) = xjp_p(:) + XB_p
        !!Write(6,"(F12.6,2x,F12.6,2x,F12.6,2x,F12.6,2x,I2,3x,I2)")  x1_p(1),x1_p(2), xj_p(1), xj_p(2), N1,N2
        !! -->  i + del_1*a_1 + del_2* a_2  =  i' + N1*L1_p + N2*L2_p  with i' in the set of lattice points.   

        
        ! The hopping increment.
        del_p  =  xj_p - xi_p

        !Twist
        pi = acos(-1.d0)
        A_p(:)  =   Flux_1 * Xnorm(Latt%a1_p) * latt%bZ1_p(:)  /  Xnorm(Latt%L1_p) + &
             &      Flux_2 * Xnorm(Latt%a2_p) * latt%bZ2_p(:)  /  Xnorm(Latt%L2_p) 

        if (Bulk) then
           !Twist in bulk
           Z_hop = Z_hop * exp(cmplx(0.d0,Iscalar(A_p,del_p),Kind(0.d0)))
        else
           !Twist as boundary
           Z_hop = Z_hop * exp(cmplx(0.d0,Iscalar(A_p,XB_p ),Kind(0.d0)))
        endif
        
        !Orbital magnetic field (Landau gauge)
        Zero =  1.0E-8
        V  =  abs(Latt%L1_p(1) * Latt%L2_p(2)  -  Latt%L1_p(2) * Latt%L2_p(1) ) 
        If ( V > Zero )  then
           B = real(N_Phi,kind(0.d0))/V
           Z_hop = Z_hop*exp(cmplx(0.d0, -2.d0*pi* B * del_p(1) *  ( xj_p(2) + xi_p(2) )/2.d0,kind(0.d0) ) )
           ! Boundary
           x_p   =  Real(N2,Kind(0.d0))*Latt%L2_p
           x1_p  =  Xjp_p + Real(N1,Kind(0.d0))*Latt%L1_p
           Z_hop =  Z_hop  *  exp(cmplx( 0.d0, -Chi(x_p, x1_p,B,pi),kind(0.d0))) 
           x_p   =  Real(N1,Kind(0.d0))*Latt%L1_p
           x1_p  =  Xjp_p 
           Z_hop =  Z_hop  *  exp(cmplx( 0.d0, -Chi(x_p, x1_p,B,pi),kind(0.d0))) 
        endif
        
        Generic_hopping =  Z_hop
        
      end function GENERIC_HOPPING

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>
!> @brief 
!> Periodic boundary conditions for Landau gauge: c_{i+L} = e{-i Chi(L,i)} c_{i}
!>    
!--------------------------------------------------------------------
      Real (Kind=kind(0.d0)) function Chi(L_p,X_p,B,pi)
        Implicit none

        Real (Kind=Kind(0.d0)), Intent(In) :: L_p(2), X_p(2), B, pi 

        Chi =  - 2.d0 * pi *B * L_p(2) * X_p(1)
      end function Chi
      
    end Module Predefined_Hoppings
