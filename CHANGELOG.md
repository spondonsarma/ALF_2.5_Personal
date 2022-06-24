# Log of backward compatibility changes and critical bugs

## 2022-06-13 Work-around for (likely) preprocessor bug

Author : J.Schwab <br>
Merge request !139

## 2022-03-21 Reset fields when no update is proposed in Global_Updates

Author : A. Goetz <br>
Merge request !136

## 2022-01-31 Write parameters to HDF5 file

Author : J.Schwab <br>
Merge request !117

### Breaking changes
1) Parameters to be formulated in format for parsing as described in Sec. 5.6 of documentation.
   Strictly speaking, it's not necessary to do that, but it simplifies the Hamiltonian,
   since the subroutine for reading parameters and writing parameters to HDF5 will be written automatically.
2) With HDF5: Add typebound procedure `write_parameters_hdf5` to Hamiltonian.

## 2021-12-08 Solves projector code runtime error

Author :  F. Parisen Toldin <br>
Merge request !129

## 2021-11-21 Implementing HDF5

Author : J.Schwab <br>
Merge request !120

### Breaking changes
1) In script configure.sh: The argument DEVEL/DEVELOPMENT is no longer a MACHINE name, but an optional switch

### Optional changes
1) Added option for compiling with HDF5 by handing argument HDF5 to configure.sh

## 2021-11-21  Automatic computation of Hopping_Matrix_Type%Multiplicity

Author : F. Parisen Toldin <br>
Merge request !116

### Breaking changes
1) Hopping_Matrix_Type%Multiplicity is now a private member, automatically initialized


## 2021-11-21  Test the checkerboard decomposition

Author : F. Parisen Toldin <br>
Merge request !124


## 2021-03-22  Implementing Submodule Hamiltonians / All hamiltonians in one binary

Author : J. Schwab <br>
Merge request !107

### Breaking changes
1) **In Hamiltonians** You will have to adapt your Hamiltonians to the Submodule structure
2) You will have to add your Hamiltonian name to the **Hamiltonians.list** in the Prog directory


## 2020-11-16   Implementing  Langevin 

Author : F.F. Assaad <br>
Merge request !91 

### Breaking changes
1) **In Hamiltonians** 

a) Mc\_step\_weight  parameter in ObserT and Obser routines <br>
b) Add 
`Subroutine Ham_Langevin_HMC_S0(Forces_0)`  <br>
Returns Bosonic forces

### Optional changes
1) **Parameters    VAR_Hubbard**

Continuous = .F.  ! Uses (T: continuous; F: discrete) HS transformation

2) **Parameters  VAR_QMC**

a) Langevin = .F.    ! Langevin update <br>
b) Delta\_t\_Langevin\_HMC = 0.01 ! Default time step for Langevin and HMC updates <br>
c) Max\_Force            = 1.5  ! Max Force for  Langevin <br>
d) HMC     = .F.   ! HMC update <br>
e) Leapfrog_steps = 0 !  Number of leapfrog steps



## 2020-09-25   Embedding lattice information in observables 

Author :  J. Schwab <br>
Merge request !66 

### Breaking changes
**In Hamiltonians** 

Calls to `Obser_Latt_make` should be adjusted to the subroutine's new interface
