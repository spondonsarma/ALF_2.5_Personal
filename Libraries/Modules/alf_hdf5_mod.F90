#if defined(HDF5)
     Module alf_hdf5
     
       interface write_attribute
         MODULE PROCEDURE write_attribute_double, write_attribute_int, write_attribute_string, write_attribute_logical
       end interface write_attribute
       interface read_attribute
         MODULE PROCEDURE read_attribute_double, read_attribute_int, read_attribute_string, read_attribute_logical
       end interface read_attribute
       interface test_attribute
         MODULE PROCEDURE test_attribute_double, test_attribute_int, test_attribute_string, test_attribute_logical
       end interface test_attribute
  
       contains
  
         Subroutine init_dset(file_id, dsetname, dims, is_complex, chunklen)
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This subroutine creates a new dataset in an opened HDF5 file for the
!> purpsose of filling if with Monte-Carlo bins.
!
!> @param [in] file_id  Idendifier of the opened HDF5 file
!> @param [in] dsetname Name of the new dataset
!> @param [in] dims     Shape of one bin. Whith size(dims) = size(bin)+1
!>                      and dims(size(dims)) = 0
!-------------------------------------------------------------------
           Use hdf5
           Implicit none
           
           INTEGER(HID_T),     intent(in) :: file_id
           Character (len=64), intent(in) :: dsetname
           INTEGER(HSIZE_T),   intent(in) :: dims(:)
           logical,            intent(in) :: is_complex
           INTEGER(HSIZE_T),   intent(in), optional :: chunklen
           
           INTEGER                       :: rank, hdferr
           INTEGER(HSIZE_T), allocatable :: dimsc(:), maxdims(:)
           INTEGER(HID_T)                :: dset_id, dataspace, crp_list
           
           !CALL h5open_f(hdferr)
           
           !Define size of dataset and of chunks
           rank = size(dims)
           allocate( dimsc(rank), maxdims(rank) )
           dimsc         = dims
           dimsc(rank)   = 1
           if ( present(chunklen) ) dimsc(rank) = chunklen
           maxdims       = dims
           maxdims(rank) = H5S_UNLIMITED_F
           
           !Check for dims(rank) = 0
           if (dims(rank) /= 0) then
             write(*,*) 'Error in init_dset: dims(rank) /= 0'
             stop
           endif
           
           !Create Dataspace
           CALL h5screate_simple_f(rank, dims, dataspace, hdferr, maxdims)
           
           !Modify dataset creation properties, i.e. enable chunking
           CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, hdferr)
           CALL h5pset_chunk_f(crp_list, rank, dimsc, hdferr)
#ifdef HDF5_ZLIB
           ! Set ZLIB / DEFLATE Compression using compression level HDF5_ZLIB
           CALL h5pset_deflate_f(crp_list, HDF5_ZLIB, hdferr)
#endif
           
           !Create a dataset using cparms creation properties.
           CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace, &
                           dset_id, hdferr, crp_list )
           
           CALL write_attribute_logical(dset_id, 'is_complex', is_complex, hdferr)
           
           !Close objects
           CALL h5sclose_f(dataspace, hdferr)
           CALL h5pclose_f(crp_list,  hdferr)
           CALL h5dclose_f(dset_id,   hdferr)
           deallocate( dimsc, maxdims )
           
         end Subroutine init_dset

!--------------------------------------------------------------------
         
         Subroutine append_dat(file_id, dsetname, dat_ptr, Nbins_in)
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This subroutine appends one bin to an existing dataset of an opened HDF5 file.
!
!> @param [in] file_id  Idendifier of the opened HDF5 file
!> @param [in] dsetname Name of the dataset
!> @param [in] data_ptr C-pointer to the first element of the data to write. 
!>                      The data should be all double precision.
!>                      The length of the data is assumed from the existing dataset.
!-------------------------------------------------------------------
           Use hdf5
           USE ISO_C_BINDING
           Implicit none
           
           INTEGER(HID_T),     intent(in) :: file_id
           Character (len=64), intent(in) :: dsetname
           TYPE(C_PTR),        intent(in) :: dat_ptr
           INTEGER, optional,  intent(in) :: Nbins_in
           
           INTEGER                       :: rank, hdferr, i, Nbins
           INTEGER(HSIZE_T)              :: mem_dims(1)
           INTEGER(HSIZE_T), allocatable :: dims(:), maxdims(:), offset(:), count(:)
           INTEGER(HID_T)                :: dset_id, dataspace, memspace
           
           !CALL h5open_f(hdferr)
           if( present(Nbins_in) ) then
             Nbins = Nbins_in
           else
             Nbins = 1
           endif
           
           !Open the  dataset.
           CALL h5dopen_f(file_id, dsetname, dset_id, hdferr)
           
           !Get dataset's dataspace handle.
           CALL h5dget_space_f(dset_id, dataspace, hdferr)
           
           !Get dataspace's rank.
           CALL h5sget_simple_extent_ndims_f(dataspace, rank, hdferr)
           allocate( dims(rank), maxdims(rank), offset(rank), count(rank) )
           
           !Get dataspace's dimensions.
           CALL h5sget_simple_extent_dims_f(dataspace, dims, maxdims, hdferr)
           
           !Extent dataset and define hyperslab to write on
           offset(:)    = 0
           offset(rank) = dims(rank)
           count(:)     = dims(:)
           count(rank)  = Nbins
           dims(rank)   = dims(rank)+Nbins
           CALL h5dset_extent_f(dset_id, dims, hdferr)
           CALL h5sclose_f(dataspace, hdferr)
           CALL h5dget_space_f(dset_id, dataspace, hdferr)
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, hdferr)
           
           !Define memory space of data
           mem_dims = Nbins
           do i=1, rank-1
             mem_dims = mem_dims*dims(i)
           enddo
           CALL h5screate_simple_f (1, mem_dims, memspace, hdferr)
           
           !Write data
           CALL H5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr, &
                           mem_space_id = memspace, file_space_id = dataspace)
           
           !Close objects
           CALL h5sclose_f(memspace, hdferr)
           CALL h5sclose_f(dataspace, hdferr)
           CALL h5dclose_f(dset_id,   hdferr)
           
           deallocate( dims, maxdims, offset, count )
           
         end Subroutine append_dat

!--------------------------------------------------------------------
         
         Subroutine write_latt(obj_id, Latt)
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This subroutine writes the lattice in an opened HDF5 object.
!
!> @param [in] obj_id  Idendifier of the opened HDF5 object
!> @param [in] Latt    The lattice
!-------------------------------------------------------------------
            Use hdf5
            use h5lt
            Use Lattices_v3
            Implicit none
            INTEGER(HID_T), intent(in) :: obj_id
            Type (Lattice), intent(in) :: Latt 
                  
            Character (len=64) :: group_name, dset_name, attr_name
            INTEGER(HID_T)  :: group_id
            LOGICAL :: link_exists
            INTEGER            :: ierr, ndim, I, rank
            INTEGER(HSIZE_T)   :: size_dat
            INTEGER(HSIZE_T), allocatable :: dims(:)
            Real (Kind=Kind(0.d0)), allocatable :: X(:,:)
            
            group_name = "lattice"
            CALL h5lexists_f(obj_id, group_name, link_exists, ierr)
            if ( link_exists ) return
            call h5gcreate_f(obj_id, group_name, group_id, ierr)
            
            ndim = size(Latt%L1_p)
            size_dat = ndim
            dset_name = "."
            attr_name = "a1"
            call h5LTset_attribute_double_f(group_id, dset_name, attr_name, Latt%a1_p, size_dat, ierr )
            attr_name = "a2"
            call h5LTset_attribute_double_f(group_id, dset_name, attr_name, Latt%a2_p, size_dat, ierr )
            attr_name = "L1"
            call h5LTset_attribute_double_f(group_id, dset_name, attr_name, Latt%L1_p, size_dat, ierr )
            attr_name = "L2"
            call h5LTset_attribute_double_f(group_id, dset_name, attr_name, Latt%L2_p, size_dat, ierr )
            !attr_name = "b1"
            !call h5LTset_attribute_double_f(group_id, dset_name, attr_name, Latt%b1_p, size_dat, ierr )
            !attr_name = "b2"
            !call h5LTset_attribute_double_f(group_id, dset_name, attr_name, Latt%b2_p, size_dat, ierr )
            !attr_name = "b1_perp"
            !call h5LTset_attribute_double_f(group_id, dset_name, attr_name, Latt%b1_perp_p, size_dat, ierr )
            !attr_name = "b2_perp"
            !call h5LTset_attribute_double_f(group_id, dset_name, attr_name, Latt%b2_perp_p, size_dat, ierr )
            !attr_name = "BZ1"
            !call h5LTset_attribute_double_f(group_id, dset_name, attr_name, Latt%BZ1_p, size_dat, ierr )
            !attr_name = "BZ2"
            !call h5LTset_attribute_double_f(group_id, dset_name, attr_name, Latt%BZ2_p, size_dat, ierr )
            
!             dset_name = "list"
!             rank = 2
!             allocate( dims(rank) )
!             dims = shape(Latt%list)
!             call h5LTmake_dataset_int_f( group_id, dset_name, rank, dims, Latt%list, ierr)
!             deallocate( dims )
!             
!             dset_name = "invlist"
!             rank = 2
!             allocate( dims(rank) )
!             dims = shape(Latt%invlist)
!             call h5LTmake_dataset_int_f( group_id, dset_name, rank, dims, Latt%invlist, ierr)
!             deallocate( dims )
!             
!             dset_name = "nnlist"
!             rank = 3
!             allocate( dims(rank) )
!             dims = shape(Latt%nnlist)
!             call h5LTmake_dataset_int_f( group_id, dset_name, rank, dims, Latt%nnlist, ierr)
!             deallocate( dims )
!             
!             dset_name = "listk"
!             rank = 2
!             allocate( dims(rank) )
!             dims = shape(Latt%listk)
!             call h5LTmake_dataset_int_f( group_id, dset_name, rank, dims, Latt%listk, ierr)
!             deallocate( dims )
!             
!             dset_name = "invlistk"
!             rank = 2
!             allocate( dims(rank) )
!             dims = shape(Latt%invlistk)
!             call h5LTmake_dataset_int_f( group_id, dset_name, rank, dims, Latt%invlistk, ierr)
!             deallocate( dims )
!             
!             dset_name = "nnlistk"
!             rank = 3
!             allocate( dims(rank) )
!             dims = shape(Latt%nnlistk)
!             call h5LTmake_dataset_int_f( group_id, dset_name, rank, dims, Latt%nnlistk, ierr)
!             deallocate( dims )
!             
!             dset_name = "imj"
!             rank = 2
!             allocate( dims(rank) )
!             dims = shape(Latt%imj)
!             call h5LTmake_dataset_int_f( group_id, dset_name, rank, dims, Latt%imj, ierr)
!             deallocate( dims )
            
            call h5gclose_f(group_id, ierr)
           
         end Subroutine write_latt

!--------------------------------------------------------------------

!--------------------------------------------------------------------

         Subroutine write_attribute_double(obj_id, attr_name, attr_value, ierr)
           Use hdf5
           Implicit none
           INTEGER(HID_T),   INTENT(IN) :: obj_id
           CHARACTER(LEN=*), INTENT(IN) :: attr_name
           real(Kind=Kind(0.d0)), INTENT(IN) :: attr_value
           INTEGER,   INTENT(OUT) :: ierr
           
           INTEGER(HID_T) :: space_id, attr_id
           INTEGER(HSIZE_T), parameter :: dims(1) = 1
           
           CALL h5screate_f (H5S_SCALAR_F, space_id, ierr)
           call h5acreate_f (obj_id, attr_name, H5T_NATIVE_DOUBLE, space_id, attr_id, ierr)
           call h5awrite_f  (attr_id, H5T_NATIVE_DOUBLE, attr_value, dims, ierr)  
           call h5aclose_f  (attr_id, ierr)
           call h5sclose_f  (space_id, ierr)
         end Subroutine write_attribute_double
        
!--------------------------------------------------------------------

         Subroutine write_attribute_int(obj_id, attr_name, attr_value, ierr)
           Use hdf5
           Implicit none
           INTEGER(HID_T),   INTENT(IN) :: obj_id
           CHARACTER(LEN=*), INTENT(IN) :: attr_name
           INTEGER,          INTENT(IN) :: attr_value
           INTEGER,   INTENT(OUT) :: ierr
           
           INTEGER(HID_T) :: space_id, attr_id
           INTEGER(HSIZE_T), parameter :: dims(1) = 1
           
           CALL h5screate_f (H5S_SCALAR_F, space_id, ierr)
           call h5acreate_f (obj_id, attr_name, H5T_NATIVE_INTEGER, space_id, attr_id, ierr)
           call h5awrite_f  (attr_id, H5T_NATIVE_INTEGER, attr_value, dims, ierr)  
           call h5aclose_f  (attr_id, ierr)
           call h5sclose_f  (space_id, ierr)
         end Subroutine write_attribute_int
        
!--------------------------------------------------------------------

         Subroutine write_attribute_string(obj_id, attr_name, attr_value, ierr)
           Use hdf5
           Implicit none
           INTEGER(HID_T),   INTENT(IN) :: obj_id
           CHARACTER(LEN=*), INTENT(IN) :: attr_name
           CHARACTER(LEN=*), INTENT(IN) :: attr_value
           INTEGER,   INTENT(OUT) :: ierr
           
           INTEGER(HID_T)     :: space_id, attr_id
           INTEGER, parameter :: rank = 1
           INTEGER(HSIZE_T)   :: dims(1)
           dims(1) = len(trim(attr_value))
           
           CALL h5screate_simple_f(rank, dims, space_id, ierr)
           call h5acreate_f (obj_id, attr_name, H5T_NATIVE_CHARACTER, space_id, attr_id, ierr)
           call h5awrite_f  (attr_id, H5T_NATIVE_CHARACTER, attr_value, dims, ierr)  
           call h5aclose_f  (attr_id, ierr)
           call h5sclose_f  (space_id, ierr)
         end Subroutine write_attribute_string
        
!--------------------------------------------------------------------

         Subroutine write_attribute_logical(obj_id, attr_name, attr_value, ierr)
           Use hdf5
           Implicit none
           INTEGER(HID_T),   INTENT(IN) :: obj_id
           CHARACTER(LEN=*), INTENT(IN) :: attr_name
           LOGICAL,          INTENT(IN) :: attr_value
           INTEGER,   INTENT(OUT) :: ierr
           
           INTEGER        :: attr_value2
           INTEGER(HID_T) :: space_id, attr_id
           INTEGER(HSIZE_T), parameter :: dims(1) = 1
           
           attr_value2 = 0
           if ( attr_value ) attr_value2 = 1
           
           CALL h5screate_f (H5S_SCALAR_F, space_id, ierr)
           call h5acreate_f (obj_id, attr_name, H5T_NATIVE_INTEGER, space_id, attr_id, ierr)
           call h5awrite_f  (attr_id, H5T_NATIVE_INTEGER, attr_value2, dims, ierr)  
           call h5aclose_f  (attr_id, ierr)
           call h5sclose_f  (space_id, ierr)
         end Subroutine write_attribute_logical
        
!--------------------------------------------------------------------

         Subroutine read_attribute_double(obj_id, attr_name, attr_value, ierr)
           Use hdf5
           Implicit none
           INTEGER(HID_T),        INTENT(IN)  :: obj_id
           CHARACTER(LEN=*),      INTENT(IN)  :: attr_name
           real(Kind=Kind(0.d0)), INTENT(OUT) :: attr_value
           INTEGER,               INTENT(OUT) :: ierr
           
           INTEGER(HID_T)              :: attr_id
           INTEGER(HSIZE_T), parameter :: dims(1) = 1
           
           call h5aopen_f  (obj_id, attr_name, attr_id, ierr)
           call h5aread_f  (attr_id, H5T_NATIVE_DOUBLE, attr_value, dims, ierr)  
           call h5aclose_f (attr_id, ierr)
         end Subroutine read_attribute_double
        
!--------------------------------------------------------------------

         Subroutine read_attribute_int(obj_id, attr_name, attr_value, ierr)
           Use hdf5
           Implicit none
           INTEGER(HID_T),   INTENT(IN)  :: obj_id
           CHARACTER(LEN=*), INTENT(IN)  :: attr_name
           INTEGER,          INTENT(OUT) :: attr_value
           INTEGER,          INTENT(OUT) :: ierr
           
           INTEGER(HID_T)              :: attr_id
           INTEGER(HSIZE_T), parameter :: dims(1) = 1
           
           call h5aopen_f  (obj_id, attr_name, attr_id, ierr)
           call h5aread_f  (attr_id, H5T_NATIVE_INTEGER, attr_value, dims, ierr)  
           call h5aclose_f (attr_id, ierr)
         end Subroutine read_attribute_int
        
!--------------------------------------------------------------------

         Subroutine read_attribute_string(obj_id, attr_name, attr_value, ierr)
           Use hdf5
           Implicit none
           INTEGER(HID_T),   INTENT(IN)  :: obj_id
           CHARACTER(LEN=*), INTENT(IN)  :: attr_name
           CHARACTER(LEN=*), INTENT(OUT) :: attr_value
           INTEGER,          INTENT(OUT) :: ierr
           
           INTEGER(HID_T)     :: attr_id
           INTEGER(HSIZE_T)   :: dims(1)
           
           call h5aopen_f  (obj_id, attr_name, attr_id, ierr)
           call h5aget_storage_size_f(attr_id, dims(1), ierr)
           attr_value = ''
           call h5aread_f  (attr_id, H5T_NATIVE_CHARACTER, attr_value, dims, ierr)
           call h5aclose_f  (attr_id, ierr)
         end Subroutine read_attribute_string
        
!--------------------------------------------------------------------

         Subroutine read_attribute_logical(obj_id, attr_name, attr_value, ierr)
           Use hdf5
           Implicit none
           INTEGER(HID_T),   INTENT(IN)  :: obj_id
           CHARACTER(LEN=*), INTENT(IN)  :: attr_name
           LOGICAL,          INTENT(OUT) :: attr_value
           INTEGER,          INTENT(OUT) :: ierr
           
           INTEGER                     :: attr_value2
           INTEGER(HID_T)              :: attr_id
           INTEGER(HSIZE_T), parameter :: dims(1) = 1
           
           call h5aopen_f  (obj_id, attr_name, attr_id, ierr)
           call h5aread_f  (attr_id, H5T_NATIVE_INTEGER, attr_value2, dims, ierr)  
           call h5aclose_f (attr_id, ierr)
           
           if ( attr_value2 == 0 ) then
             attr_value = .false.
           elseif ( attr_value2 == 1 ) then
             attr_value = .true.
           else
             write(*,*) "Error in read_attribute_logical: attr_value2 is neither 0 or 1, but", attr_value2
             stop
           endif
         end Subroutine read_attribute_logical
        
!--------------------------------------------------------------------

         Subroutine test_attribute_double(obj_id, attr_name, attr_value, ierr)
           Use hdf5
           Implicit none
           INTEGER(HID_T),   INTENT(IN) :: obj_id
           CHARACTER(LEN=*), INTENT(IN) :: attr_name
           real(Kind=Kind(0.d0)), INTENT(IN) :: attr_value
           INTEGER,   INTENT(OUT) :: ierr
           
           LOGICAL :: attr_exists
           real(Kind=Kind(0.d0)), parameter :: ZERO = 10D-8
           real(Kind=Kind(0.d0)) :: test_double, diff
           
           call h5aexists_f(obj_id, attr_name, attr_exists, ierr)
           
           if ( .not. attr_exists ) then
             call write_attribute_double(obj_id, attr_name, attr_value, ierr)
           else
             call read_attribute_double(obj_id, attr_name, test_double, ierr)
             diff = abs(attr_value - test_double)
             if (diff > ZERO) then
               write(*,*) 'Error in test_attribute_double:', attr_name, ' = ', attr_value, '/=', test_double
               stop
             endif
           endif
         end Subroutine test_attribute_double
        
!--------------------------------------------------------------------

         Subroutine test_attribute_int(obj_id, attr_name, attr_value, ierr)
           Use hdf5
           Implicit none
           INTEGER(HID_T),   INTENT(IN) :: obj_id
           CHARACTER(LEN=*), INTENT(IN) :: attr_name
           INTEGER,          INTENT(IN) :: attr_value
           INTEGER,   INTENT(OUT) :: ierr
           
           LOGICAL :: attr_exists
           INTEGER :: test_int
           
           call h5aexists_f(obj_id, attr_name, attr_exists, ierr)
           
           if ( .not. attr_exists ) then
             call write_attribute_int(obj_id, attr_name, attr_value, ierr)
           else
             call read_attribute_int(obj_id, attr_name, test_int, ierr)
             if (attr_value /= test_int) then
               write(*,*) 'Error in test_attribute_int:', attr_name, ' = ', attr_value, '/=', test_int
               stop
             endif
           endif
         end Subroutine test_attribute_int
        
!--------------------------------------------------------------------

         Subroutine test_attribute_string(obj_id, attr_name, attr_value, ierr)
           Use hdf5
           Implicit none
           INTEGER(HID_T),   INTENT(IN) :: obj_id
           CHARACTER(LEN=*), INTENT(IN) :: attr_name
           CHARACTER(LEN=*), INTENT(IN) :: attr_value
           INTEGER,   INTENT(OUT) :: ierr
           
           LOGICAL :: attr_exists
           CHARACTER(LEN=64) :: test_string
           
           call h5aexists_f(obj_id, attr_name, attr_exists, ierr)
           
           if ( .not. attr_exists ) then
             call write_attribute_string(obj_id, attr_name, attr_value, ierr)
           else
             call read_attribute_string(obj_id, attr_name, test_string, ierr)
             if (trim(attr_value) /= trim(test_string)) then
               write(*,*) 'Error in test_attribute_string:', attr_name, ' = ', attr_value, '/=', test_string
               stop
             endif
           endif
         end Subroutine test_attribute_string
        
!--------------------------------------------------------------------

         Subroutine test_attribute_logical(obj_id, attr_name, attr_value, ierr)
           Use hdf5
           Implicit none
           INTEGER(HID_T),   INTENT(IN) :: obj_id
           CHARACTER(LEN=*), INTENT(IN) :: attr_name
           LOGICAL,          INTENT(IN) :: attr_value
           INTEGER,   INTENT(OUT) :: ierr
           
           LOGICAL :: attr_exists
           LOGICAL :: test_bool
           
           call h5aexists_f(obj_id, attr_name, attr_exists, ierr)
           
           if ( .not. attr_exists ) then
             call write_attribute_logical(obj_id, attr_name, attr_value, ierr)
           else
             call read_attribute_logical(obj_id, attr_name, test_bool, ierr)
             if (attr_value .neqv. test_bool) then
               write(*,*) 'Error in test_attribute_logical:', attr_name, ' = ', attr_value, '/=', test_bool
               stop
             endif
           endif
         end Subroutine test_attribute_logical
         
     end Module alf_hdf5
#endif
