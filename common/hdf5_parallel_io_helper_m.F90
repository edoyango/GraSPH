module hdf5_parallel_io_helper_m
   use param, only: f, inputf, outputf
   use hdf5
#ifdef PARALLEL
   use mpi ! would prefer to use mpi_f08, but hdf5 doesn't play nicely with it
#endif
   private

   integer(HID_T):: plist_id, compress_plist_id, dspace_id, dset_id, local_dspace_id
   integer:: ierr

   interface hdf5_parallel_read
      module procedure hdf5_parallel_read_flt_r1, hdf5_parallel_read_flt_r2, hdf5_parallel_read_int_r1
   end interface hdf5_parallel_read

   interface hdf5_parallel_write
      module procedure hdf5_parallel_write_flt_r1, hdf5_parallel_write_flt_r2, hdf5_parallel_write_int_r1
   end interface hdf5_parallel_write

   public:: hdf5_fileopen_read, hdf5_fileopen_write, hdf5_attribute_write, hdf5_attribute_read, hdf5_parallel_read, &
      hdf5_parallel_write

contains

   !===============================================================================================================================
   subroutine hdf5_fileopen_write(fname, fid)

      implicit none
      character(len=*):: fname
      integer(HID_T), intent(out):: fid

      ! creating file access property list
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
#ifdef PARALLEL
      call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
#endif
      call h5fcreate_f(fname, H5F_ACC_TRUNC_F, fid, ierr, access_prp=plist_id)
      call h5pclose_f(plist_id, ierr)

   end subroutine hdf5_fileopen_write

   !===============================================================================================================================
   subroutine hdf5_fileopen_read(fname, fid)

      implicit none
      character(len=*):: fname
      integer(HID_T), intent(out):: fid

      ! creating file access property list
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
#ifdef PARALLEL
      call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
#endif
      call h5fopen_f(fname, H5F_ACC_RDONLY_F, fid, ierr, access_prp=plist_id)
      call h5pclose_f(plist_id, ierr)

   end subroutine hdf5_fileopen_read

   !====================================================================================================================
   subroutine hdf5_attribute_write(gid, n)

      implicit none
      integer(HID_T), intent(in):: gid
      integer, intent(in):: n
      integer(HID_T):: aspace_id, attr_id
      integer(SIZE_T):: dims(1) = 0
      integer:: ierr

      call h5screate_f(H5S_SCALAR_F, aspace_id, ierr)
      call h5acreate_f(gid, 'n', H5T_NATIVE_INTEGER, aspace_id, attr_id, ierr)
      call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, n, dims, ierr)
      call h5aclose_f(attr_id, ierr)
      call h5sclose_f(aspace_id, ierr)

   end subroutine hdf5_attribute_write 

   !====================================================================================================================
   subroutine hdf5_attribute_read(gid, n)

      implicit none
      integer(HID_T), intent(in):: gid
      integer, intent(out):: n
      integer(HID_T):: aspace_id, attr_id
      integer(size_t):: dims(1) = 0
      integer:: ierr

      call h5aopen_by_name_f(gid, '.', 'n', attr_id, ierr)
      call h5aget_space_f(attr_id, aspace_id, ierr)
      call h5aread_f(attr_id, H5T_NATIVE_INTEGER, n, dims, ierr)
      call h5sclose_f(aspace_id, ierr)
      call h5aclose_f(attr_id, ierr)

   end subroutine hdf5_attribute_read

! #ifdef PARALLEL
   !===============================================================================================================================
   subroutine hdf5_parallel_write_flt_r1(gid, dset_name, displ, global_dims, ddata)

      implicit none
      integer(HID_T), intent(in):: gid
      character(len=*), intent(in):: dset_name
      integer(HSIZE_T), intent(in):: global_dims ! shape of entire dataset being written
      integer(HSSIZE_T), intent(in):: displ ! displacement of process
      real(outputf), intent(in):: ddata(:) ! data local to MPI Process to write
      integer(HSIZE_T):: local_dims(1), global_dims_copy(1), chunk_size(1)
      integer(HSSIZE_T):: displ_copy(1)
      integer, parameter:: rank = 1

      local_dims = size(ddata)
      global_dims_copy(1) = global_dims
      displ_copy(1) = displ

      ! Creating dataspace for entire dataset
      call h5screate_simple_f(rank, global_dims_copy, dspace_id, ierr)

      chunk_size = global_dims_copy  ! need a better way to choose a value
      call h5pcreate_f(H5P_DATASET_CREATE_F, compress_plist_id, ierr)
      call h5pset_chunk_f(compress_plist_id, 1, chunk_size, ierr)
      call h5pset_deflate_f(compress_plist_id, 6, ierr)

      ! Creating dataset for entire dataset
      select case (outputf)
      case (kind(1.))
         call h5dcreate_f(gid, dset_name, H5T_NATIVE_REAL, dspace_id, dset_id, ierr, dcpl_id=compress_plist_id)
      case default
         call h5dcreate_f(gid, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr, dcpl_id=compress_plist_id)
      end select
      call h5sclose_f(dspace_id, ierr)

      ! Creating dataspace in dataset for each process to write to
      call h5screate_simple_f(rank, local_dims, local_dspace_id, ierr)
      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, displ_copy, local_dims, ierr)
#ifdef PARALLEL
      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
#else
      plist_id = H5P_DEFAULT_F
#endif
      ! Write the dataset collectively.
      select case (outputf)
      case (kind(1.))
         call h5dwrite_f(dset_id, H5T_NATIVE_REAL, ddata, global_dims_copy, ierr, file_space_id=dspace_id, &
                         mem_space_id=local_dspace_id, xfer_prp=plist_id)
      case default
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ddata, global_dims_copy, ierr, file_space_id=dspace_id, &
                         mem_space_id=local_dspace_id, xfer_prp=plist_id)
      end select

      ! closing all dataspaces,datasets,property lists.
      call h5sclose_f(dspace_id, ierr)
      call h5sclose_f(local_dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
      call h5pclose_f(plist_id, ierr)
      call h5pclose_f(compress_plist_id, ierr)

   end subroutine hdf5_parallel_write_flt_r1

   !===============================================================================================================================
   subroutine hdf5_parallel_read_flt_r1(gid, dset_name, displ, global_dims, ddata)

      implicit none
      integer(HID_T), intent(in):: gid
      character(len=*), intent(in):: dset_name
      integer(HSIZE_T), intent(in):: global_dims ! shape of entire dataset being written
      integer(HSSIZE_T), intent(in):: displ ! displacement of process
      real(f), intent(out):: ddata(:) ! data local to MPI Process to write
      integer(HSIZE_T):: local_dims(1), global_dims_copy(1)
      integer(HSSIZE_T):: displ_copy(1)
      integer, parameter:: rank = 1

      local_dims = size(ddata)
      global_dims_copy(1) = global_dims
      displ_copy(1) = displ

      ! Creating dataspace for entire dataset
      call h5screate_simple_f(rank, global_dims_copy, dspace_id, ierr)

      ! Creating dataset for entire dataset
      call h5dopen_f(gid, dset_name, dset_id, ierr)
      call h5sclose_f(dspace_id, ierr)

      ! Creating dataspace in dataset for each process to write to
      call h5screate_simple_f(rank, local_dims, local_dspace_id, ierr)
      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, displ_copy, local_dims, ierr)

#ifdef PARALLEL
      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
#else
      plist_id = H5P_DEFAULT_F
#endif
      ! Write the dataset collectively.
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, ddata, global_dims_copy, ierr, file_space_id=dspace_id, &
         mem_space_id=local_dspace_id, xfer_prp=plist_id)

      ! closing all dataspaces,datasets,property lists.
      call h5sclose_f(dspace_id, ierr)
      call h5sclose_f(local_dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
      call h5pclose_f(plist_id, ierr)

   end subroutine hdf5_parallel_read_flt_r1

   !===============================================================================================================================
   subroutine hdf5_parallel_write_flt_r2(gid, dset_name, displ, global_dims, ddata)

      implicit none
      integer(HID_T), intent(in):: gid
      character(len=*), intent(in):: dset_name
      integer(HSIZE_T), intent(in):: global_dims(2) ! shape of entire dataset being written
      integer(HSSIZE_T), intent(in):: displ(2) ! displacement of process
      real(outputf), intent(in):: ddata(:, :) ! data local to MPI Process to write
      integer(HSIZE_T):: local_dims(2), chunk_size(2)
      integer, parameter:: rank = 2

      local_dims = shape(ddata)

      ! Creating dataspace for entire dataset
      call h5screate_simple_f(rank, global_dims, dspace_id, ierr)

      chunk_size(1) = global_dims(1)  ! need a better way to choose a value
      chunk_size(2) = global_dims(2)
      call h5pcreate_f(H5P_DATASET_CREATE_F, compress_plist_id, ierr)
      call h5pset_chunk_f(compress_plist_id, 2, chunk_size, ierr)
      call h5pset_deflate_f(compress_plist_id, 6, ierr)

      ! Creating dataset for entire dataset
      select case (outputf)
      case (kind(1.))
         call h5dcreate_f(gid, dset_name, H5T_NATIVE_REAL, dspace_id, dset_id, ierr, dcpl_id=compress_plist_id)
      case default
         call h5dcreate_f(gid, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr, dcpl_id=compress_plist_id)
      end select
      call h5sclose_f(dspace_id, ierr)

      ! Creating dataspace in dataset for each process to write to
      call h5screate_simple_f(rank, local_dims, local_dspace_id, ierr)
      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, displ, local_dims, ierr)
#ifdef PARALLEL
      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
#else
      plist_id = H5P_DEFAULT_F
#endif
      ! Write the dataset collectively.
      select case (outputf)
      case (kind(1.))
         call h5dwrite_f(dset_id, H5T_NATIVE_REAL, ddata, global_dims, ierr, file_space_id=dspace_id, &
                         mem_space_id=local_dspace_id, xfer_prp=plist_id)
      case default
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ddata, global_dims, ierr, file_space_id=dspace_id, &
                         mem_space_id=local_dspace_id, xfer_prp=plist_id)
      end select

      ! closing all dataspaces,datasets,property lists.
      call h5sclose_f(dspace_id, ierr)
      call h5sclose_f(local_dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
      call h5pclose_f(plist_id, ierr)
      call h5pclose_f(compress_plist_id, ierr)

   end subroutine hdf5_parallel_write_flt_r2

   !===============================================================================================================================
   subroutine hdf5_parallel_read_flt_r2(gid, dset_name, displ, global_dims, ddata)

      implicit none
      integer(HID_T), intent(in):: gid
      character(len=*), intent(in):: dset_name
      integer(HSIZE_T), intent(in):: global_dims(2) ! shape of entire dataset being written
      integer(HSSIZE_T), intent(in):: displ(2) ! displacement of process
      real(f), intent(out):: ddata(:, :) ! data local to MPI Process to write
      integer(HSIZE_T):: local_dims(2)
      integer, parameter:: rank = 2

      local_dims = shape(ddata)

      ! Creating dataspace for entire dataset
      call h5screate_simple_f(rank, global_dims, dspace_id, ierr)

      ! Creating dataset for entire dataset
      call h5dopen_f(gid, dset_name, dset_id, ierr)
      call h5sclose_f(dspace_id, ierr)

      ! Creating dataspace in dataset for each process to write to
      call h5screate_simple_f(rank, local_dims, local_dspace_id, ierr)
      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, displ, local_dims, ierr)
#ifdef PARALLEL
      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
#else
      plist_id = H5P_DEFAULT_F
#endif
      ! Write the dataset collectively.
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, ddata, global_dims, ierr, file_space_id=dspace_id, &
         mem_space_id=local_dspace_id, xfer_prp=plist_id)

      ! closing all dataspaces,datasets,property lists.
      call h5sclose_f(dspace_id, ierr)
      call h5sclose_f(local_dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
      call h5pclose_f(plist_id, ierr)

   end subroutine hdf5_parallel_read_flt_r2

   !===============================================================================================================================
   subroutine hdf5_parallel_write_int_r1(gid, dset_name, displ, global_dims, idata)

      implicit none
      integer(HID_T), intent(in):: gid
      character(len=*), intent(in):: dset_name
      integer(HSIZE_T), intent(in):: global_dims ! shape of entire dataset being written
      integer(HSSIZE_T), intent(in):: displ ! displacement of process
      integer, intent(in):: idata(:) ! data local to MPI Process to write
      integer(HSIZE_T):: local_dims(1), global_dims_copy(1), chunk_size(1)
      integer(HSSIZE_T):: displ_copy(1) ! this exists because the hdf5 subroutine only accepts arrays, not constants
      integer, parameter:: rank = 1

      local_dims = size(idata)
      global_dims_copy(1) = global_dims
      displ_copy(1) = displ

      ! Creating dataspace for entire dataset
      call h5screate_simple_f(rank, global_dims_copy, dspace_id, ierr)

      chunk_size = global_dims_copy ! need a better way to choose a value
      call h5pcreate_f(H5P_DATASET_CREATE_F, compress_plist_id, ierr)
      call h5pset_chunk_f(compress_plist_id, 1, chunk_size, ierr)
      call h5pset_deflate_f(compress_plist_id, 6, ierr)

      ! Creating dataset for entire dataset
      call h5dcreate_f(gid, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr, dcpl_id=compress_plist_id)
      call h5sclose_f(dspace_id, ierr)

      ! Creating dataspace in dataset for each process to write to
      call h5screate_simple_f(rank, local_dims, local_dspace_id, ierr)
      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, displ_copy, local_dims, ierr)
#ifdef PARALLEL
      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
#else
      plist_id = H5P_DEFAULT_F
#endif
      ! Write the dataset collectively.
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, idata, global_dims_copy, ierr, file_space_id=dspace_id, &
                      mem_space_id=local_dspace_id, xfer_prp=plist_id)

      ! closing all dataspaces,datasets,property lists.
      call h5sclose_f(dspace_id, ierr)
      call h5sclose_f(local_dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
      call h5pclose_f(plist_id, ierr)
      call h5pclose_f(compress_plist_id, ierr)

   end subroutine hdf5_parallel_write_int_r1

   !===============================================================================================================================
   subroutine hdf5_parallel_read_int_r1(gid, dset_name, displ, global_dims, idata)

      implicit none
      integer(HID_T), intent(in):: gid
      character(len=*), intent(in):: dset_name
      integer(HSIZE_T), intent(in):: global_dims ! shape of entire dataset being written
      integer(HSSIZE_T), intent(in):: displ ! displacement of process
      integer, intent(out):: idata(:) ! data local to MPI Process to write
      integer(HSIZE_T):: local_dims(1), global_dims_copy(1)
      integer(HSSIZE_T):: displ_copy(1)
      integer, parameter:: rank = 1

      local_dims = size(idata)
      global_dims_copy(1) = global_dims
      displ_copy(1) = displ

      ! Creating dataspace for entire dataset
      call h5screate_simple_f(rank, global_dims_copy, dspace_id, ierr)

      ! Creating dataset for entire dataset
      call h5dopen_f(gid, dset_name, dset_id, ierr)
      call h5sclose_f(dspace_id, ierr)

      ! Creating dataspace in dataset for each process to write to
      call h5screate_simple_f(rank, local_dims, local_dspace_id, ierr)
      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, displ_copy, local_dims, ierr)
#ifdef PARALLEL
      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
#else
      plist_id = H5P_DEFAULT_F
#endif
      ! Write the dataset collectively.
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, idata, global_dims_copy, ierr, file_space_id=dspace_id, &
         mem_space_id=local_dspace_id, xfer_prp=plist_id)

      ! closing all dataspaces,datasets,property lists.
      call h5sclose_f(dspace_id, ierr)
      call h5sclose_f(local_dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
      call h5pclose_f(plist_id, ierr)

   end subroutine hdf5_parallel_read_int_r1
! #endif
end module hdf5_parallel_io_helper_m
