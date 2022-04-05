program check_libs

  use iso_fortran_env
  integer :: N=2, ix=1, iy=1,  LDA=2, LWMAX=1000, INFO, LWORK=-1
  real    (kind=kind(0.d0)) :: C(2), D(2), dot,   W(2), RWORK(4)
  complex (kind=kind(0.d0)) :: A(2,2), WORK(1000)
  real    (kind=kind(0.d0)), external :: ddot
  data    A /(1,0),(0,0),(0,0),(2,0)/,  C /3,1415/, D /42,0/
  external ZHEEV
  
  dot = ddot(n,C,ix,D,iy)
  call ZHEEV( 'V', 'U', N, A, LDA, W, WORK, LWORK, RWORK, INFO )
  
end program check_libs
