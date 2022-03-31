program check_libs
  use iso_fortran_env
  integer :: n,ix=1,iy=1
  real(kind=8) :: A(1),B(1)
  real(kind=8) :: z
  real(kind=8), external :: ddot
  A(1)=1; B(1)=1
  z=ddot(n,A,ix,B,iy)
  !write(*,*) "Linear algebra libraries found."
end program check_libs
