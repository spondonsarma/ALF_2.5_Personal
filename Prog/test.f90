program test
    Use DynamicMatrixArray_mod
    Use ContainerElementBase_mod
    Use OpTTypes_mod
    Use Operator_mod
    implicit none

    Type(DynamicMatrixArray) :: vec
    Type(Operator), dimension(:), allocatable :: Op_T
    Class(RealExpOpT), pointer :: reopt => null()
    Class(CmplxExpOpT), pointer :: cmplxopt => null()
    Class(ContainerElementBase), pointer :: dummy
    Complex(kind=kind(0.d0)), allocatable, dimension(:,:) :: res
    Complex(kind=kind(0.d0)) :: alpha, zero
    Integer :: i, j, k, l, nmax, ndimmax

    nmax = 5
    ndimmax = 5
    allocate (res(ndimmax, ndimmax))
    call vec%init()

    alpha = 1.0
    zero = 0.0
    ! initialize res as unit matrix
    call zlaset('A', ndimmax, ndimmax, zero, alpha, res, ndimmax)

    allocate(Op_T(2*nmax))

    do i = 1, nmax
    Call Op_make(Op_T(i), Ndimmax)
    Call Op_make(Op_T(nmax + i), Ndimmax)
    Op_T(i)%O = 0
    Op_T(nmax + i)%O = 0
    Do j = 1, ndimmax
    Op_T(i)%P(j) = j
    Op_T(nmax + i)%P(j) = j
    Enddo
    Op_T(i)%g      = 1
    Op_T(i)%type      = 2
    Op_T(i)%alpha = cmplx(0.d0,0.d0, kind(0.D0))
    Op_T(nmax + i)%g      = 1
    Op_T(nmax + i)%type      = 2
    Op_T(nmax + i)%alpha = cmplx(0.d0,0.d0, kind(0.D0))
    ! fill with some data
    do j = 1, ndimmax
    Op_T(i)%O(j,j) = j
    Op_T(nmax + i)%O(j,j) = cmplx(j, j, kind(0.D0))
    enddo

    Call Op_set(Op_T(i))
    Call Op_set(Op_T(nmax + i))

    allocate(reopt, cmplxopt)
    call reopt%init(Op_T(i))
    call cmplxopt%init(Op_T(nmax + i))
    call vec%pushback(reopt)
    call vec%pushback(cmplxopt)
    enddo


    ! execute a loop over all stored objects
    do i= 1, vec%length()
    dummy => vec%at(i) ! get object
    call dummy%rmult(res) ! polymorphic dispatch to rmult
!     do k = 1, nmax
!     write (*,*) (dble(res(k,l)), l = 1,nmax )
!     enddo
!     write (*,*) "============"
    enddo

    do i = 1, nmax
    write (*,*) (res(i,j), j = 1,nmax )
    enddo

    ! tidy up
    do i = 1, vec%length()
    dummy => vec%at(i) ! Fortran doesn't want chaining here
    deallocate(dummy)
    call Op_clear(Op_T(i), ndimmax)
    enddo
    call vec%dealloc()
    deallocate(res, Op_T)
end program
