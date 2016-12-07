Program UDVC
Use MyMats
implicit none
COMPLEX (KIND=KIND(0.D0)), dimension(:, :), allocatable :: a
COMPLEX (KIND=KIND(0.D0)), dimension(:, :), allocatable :: u, v, test
COMPLEX (KIND=KIND(0.D0)), dimension(:), allocatable :: d
INTEGER :: ncon, i, j, ndim, nr
COMPLEX (KIND=KIND(0.D0)) :: Z

ndim = 5
ncon = 1
allocate(a(ndim, ndim), u(ndim, ndim), v(ndim, ndim), d(ndim), test(ndim, ndim))
! initialize A to the Hilbert matrix
do i = 1, ndim
do j = 1, ndim
a(i, j) = 1.D0/(i+j)
enddo
enddo
call UDV_C(a, u, d, v, ncon)
! check decomposition
           DO J = 1,ndim
              DO I = 1,ndim
                 Z = 0.D0
                 DO NR = 1,ndim 
                    Z = Z + U(I,NR)*D(NR)*V(NR,J)
                 ENDDO
                 TEST(I,J) = Z
              ENDDO
           ENDDO
deallocate(a, u, d, v, test)
end program
