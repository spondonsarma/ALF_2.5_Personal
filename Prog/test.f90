program test
Use DynamicMatrixArray_mod
Use ContainerElementBase_mod
Use matTypes_mod
implicit none

Type(DynamicMatrixArray) :: vec
Type(RealMat), allocatable :: remat
Type(CmplxMat), allocatable:: complexmat
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

allocate(remat, complexmat)

do i = 1, 5
    ! create some complex dummy data
    call zlaset('A', nmax, nmax, zero, alpha, ctmp, nmax)
    do j = 1, nmax
    ctmp(j,j) = i
    enddo

    !pushback
    call complexmat%init(ctmp)
    call vec%pushback(complexmat)

    ! create some real dummy data
    call dlaset('A', nmax, nmax, zero, alpha, rtmp, nmax)
    do j = 1, nmax
    rtmp(j,j) = i+j
    enddo
    ! push_back
    call remat%init(rtmp)
    call vec%pushback(remat)
enddo
! tidy up auxiliary structures
deallocate(remat, complexmat)
deallocate(ctmp, rtmp)

! execute a loop over all stored objects
do i= 1, vec%length()
    dummy = vec%at(i) ! get object
    call dummy%rmult(res) ! polymorphic dispatch to rmult
    do k = 1, nmax
    write (*,*) (dble(res(k,l)), l = 1,nmax )
    enddo
    write (*,*) "============"
enddo

do i = 1, nmax
write (*,*) (res(i,j), j = 1,nmax )
enddo

! tidy up
do i = 1, vec%length()
dummy = vec%at(i) ! Fortran doesn't want chaining here
deallocate(dummy)
enddo
call vec%dealloc()
deallocate(res)
end program
