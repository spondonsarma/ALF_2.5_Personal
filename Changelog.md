Implementing  Langevin: 

1) **Parameters    VAR_Hubbard**

Continuous = .F.  ! Uses (T: continuous; F: discrete) HS transformation

2) **Parameters  VAR_QMC**

a) Global\_update\_scheme = "Langevin"   ! Langevin or HMC <br>
b) Delta\_t\_Langevin_HMC = 0.01 ! Default time step for Langevin and HMC updates <br>
c) Max\_Force            = 1.5  ! Max Force for  Langevin

3) **In Hamiltonians**

a) Mc\_step\_weight  parameter in ObserT and Obser <br>
b) Subroutine Ham\_Langevin\_HMC\_S0(Forces\_0)  <br>
Returns Bosonic forces