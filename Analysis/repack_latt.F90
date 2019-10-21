Program repack_latt
  use ana_mod
  Implicit none
  
  INTEGER :: i, nb, no, no1, nt
  Character (len=64) :: File_in, File_out !, groupname, obs_dsetname, bak_dsetname, sgn_dsetname
!   INTEGER(HID_T)     :: file_id, group_id
!   logical            :: link_exists
!   INTEGER(HSIZE_T), allocatable :: dims(:)
!   TYPE(C_PTR)                   :: dat_ptr
!   real(Kind=Kind(0.d0)), target :: sgn
  
  Real    (Kind=Kind(0.d0)), allocatable :: phase(:)
  Complex (Kind=Kind(0.d0)), pointer     :: bins(:,:,:,:,:), bins0(:,:)
  Complex (Kind=Kind(0.d0)), pointer     :: bin(:,:,:,:), bin0(:)
  Integer                                :: Norb, Nunit, Ntau, Nbins
  Real    (Kind=Kind(0.d0))              :: dtau, X_p(2)
  logical                                :: timedisplaced
  Type (Lattice)    :: Latt
  
  i = 1
  CALL GET_COMMAND_ARGUMENT(i, File_in)
  i = 2
  CALL GET_COMMAND_ARGUMENT(i, File_out)
  
  write(*,*) "reading from ", File_in
  !determining whether equal time or tau. sloppy
  i = len(trim(File_in)) -2
  if ( File_in(i:) == 'tau' ) then
    timedisplaced = .true.
  else
    timedisplaced = .false.
  endif
  write(*,*) timedisplaced
  
  if ( timedisplaced ) then
    call read_latt(File_in, phase, bins, bins0, Latt, dtau)
  else
    call read_latt(File_in, phase, bins, bins0, Latt)
    dtau = -1.d0
  endif
  Nunit  = size(bins, 1) 
  Ntau   = size(bins, 2)
  Norb   = size(bins, 3)
  Nbins  = size(bins, 5)
  
  write(*,*) "writing to ", File_out
  
  Open (Unit=10,File=File_out, status="unknown",  position="append")
  do nb = 1,Nbins
    Write(10,*) phase(nb), Norb, Latt%N, Ntau, dtau, &
                & size(Latt%L1_p), Latt%L1_p(:), Latt%L2_p(:), Latt%a1_p(:), Latt%a2_p(:)
    Do no = 1,Norb
      Write(10,*) bins0(no,nb)
    enddo
    do I = 1,Latt%N
      x_p = dble(Latt%listk(i,1))*Latt%b1_p + dble(Latt%listk(i,2))*Latt%b2_p  
      Write(10,*) X_p(1), X_p(2)
      Do nt = 1,Ntau
        do no = 1,Norb
          do no1 = 1,Norb
            Write(10,*) bins(I,nt,no,no1,nb)
          enddo
        enddo
      enddo
    enddo
  enddo
  close(10)
  
  
end Program
