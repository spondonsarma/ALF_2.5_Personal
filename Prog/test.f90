program test
Use DynamicMatrixArray_mod
Use ContainerElementBase_mod
Use OpTTypes_mod
implicit none

Type(DynamicMatrixArray) :: vec
Type(RealOpT), allocatable :: remat
Type(CmplxOpT), allocatable:: cmplxmat
class(ContainerElementBase), allocatable :: dummy
Complex(kind=kind(0.d0)), allocatable, dimension(:,:) :: res, ctmp
Real(kind=kind(0.d0)), allocatable, dimension(:,:) :: rtmp
Complex(kind=kind(0.d0)) :: alpha, zero
Integer :: i,j,k,l, nmax

nmax = 5
allocate (res(nmax, nmax), ctmp(nmax, nmax), rtmp(nmax, nmax))
call vec%init()

alpha = 1.0
zero = 0.0
call zlaset('A', nmax, nmax, zero, alpha, res, nmax)

allocate(remat, cmplxmat)

do i = 1, 5
    call zlaset('A', nmax, nmax, zero, alpha, ctmp, nmax)
    do j = 1, nmax
    ctmp(j,j) = j
    enddo

    call cmplxmat%init(ctmp)
    call vec%pushback(cmplxmat)

    call dlaset('A', nmax, nmax, zero, alpha, rtmp, nmax)
    do j = 1, nmax
    rtmp(j,j) = j
    enddo
    call remat%init(rtmp)
    call vec%pushback(remat)
enddo
deallocate(remat, cmplxmat)
deallocate(ctmp, rtmp)

do i= 1, 5
call vec%at(i, dummy)
call dummy%rmult(res)
do k = 1, nmax
write (*,*) (res(k,l), l = 1,nmax )
enddo
write (*,*) "============"
enddo

do i = 1, nmax
write (*,*) (res(i,j), j = 1,nmax )
enddo

! tidy up
do i = 1, vec%length()
call vec%at(i, dummy)
deallocate(dummy)
enddo
call vec%dealloc()
deallocate(res)
end program
