module hdf5_parallel_io_helper_m
#ifdef PARALLEL
   use hdf5
   use mpi ! would prefer to use mpi_f08, but hdf5 doesn't play nicely with it
   use param, only: f

   interface hdf5_parallel_write
      module procedure hdf5_parallel_write_flt_r1, hdf5_parallel_write_flt_r2, hdf5_parallel_write_int_r1
   end interface hdf5_parallel_write

   integer(HID_T), private:: plist_id, dspace_id, dset_id, local_dspace_id
   integer, private:: ierr

   public:: hdf5_parallel_write, hdf5_parallel_fileopen
   private::  hdf5_parallel_write_flt_r1, hdf5_parallel_write_flt_r2, hdf5_parallel_write_int_r1
#endif
contains
#ifdef PARALLEL
   !===============================================================================================================================
   subroutine hdf5_parallel_fileopen(fname, fid)

      implicit none
      character(len=*):: fname
      integer(HID_T), intent(out):: fid

      ! creating file access property list
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)

      call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
      call h5fcreate_f(fname, H5F_ACC_TRUNC_F, fid, ierr, access_prp=plist_id)
      call h5pclose_f(plist_id, ierr)

   end subroutine hdf5_parallel_fileopen

   !===============================================================================================================================
   subroutine hdf5_parallel_write_flt_r1(fid, dset_name, displ, global_dims, ddata)

      implicit none
      integer(HID_T), intent(in):: fid
      character(len=*), intent(in):: dset_name
      integer(HSIZE_T), intent(in):: global_dims ! shape of entire dataset being written
      integer(HSSIZE_T), intent(in):: displ ! displacement of process
      real(f), intent(in):: ddata(:) ! data local to MPI Process to write
      integer(HSIZE_T):: local_dims(1), global_dims_copy(1)
      integer(HSSIZE_T):: displ_copy(1)
      integer, parameter:: rank = 1

      local_dims = size(ddata)
      global_dims_copy(1) = global_dims
      displ_copy(1) = displ

      ! Creating dataspace for entire dataset
      call h5screate_simple_f(rank, global_dims_copy, dspace_id, ierr)

      ! Creating dataset for entire dataset
      call h5dcreate_f(fid, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
      call h5sclose_f(dspace_id, ierr)

      ! Creating dataspace in dataset for each process to write to
      call h5screate_simple_f(rank, local_dims, local_dspace_id, ierr)
      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, displ_copy, local_dims, ierr)
      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
      ! Write the dataset collectively.
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ddata, global_dims_copy, ierr, file_space_id=dspace_id, &
                      mem_space_id=local_dspace_id, xfer_prp=plist_id)

      ! closing all dataspaces,datasets,property lists.
      call h5sclose_f(dspace_id, ierr)
      call h5sclose_f(local_dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
      call h5pclose_f(plist_id, ierr)

   end subroutine hdf5_parallel_write_flt_r1

   !===============================================================================================================================
   subroutine hdf5_parallel_write_flt_r2(fid, dset_name, displ, global_dims, ddata)

      implicit none
      integer(HID_T), intent(in):: fid
      character(len=*), intent(in):: dset_name
      integer(HSIZE_T), intent(in):: global_dims(2) ! shape of entire dataset being written
      integer(HSSIZE_T), intent(in):: displ(2) ! displacement of process
      real(f), intent(in):: ddata(:, :) ! data local to MPI Process to write
      integer(HSIZE_T):: local_dims(2)
      integer, parameter:: rank = 2

      local_dims = shape(ddata)

      ! Creating dataspace for entire dataset
      call h5screate_simple_f(rank, global_dims, dspace_id, ierr)

      ! Creating dataset for entire dataset
      call h5dcreate_f(fid, dset_name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
      call h5sclose_f(dspace_id, ierr)

      ! Creating dataspace in dataset for each process to write to
      call h5screate_simple_f(rank, local_dims, local_dspace_id, ierr)
      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, displ, local_dims, ierr)
      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
      ! Write the dataset collectively.
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ddata, global_dims, ierr, file_space_id=dspace_id, mem_space_id=local_dspace_id, &
                      xfer_prp=plist_id)

      ! closing all dataspaces,datasets,property lists.
      call h5sclose_f(dspace_id, ierr)
      call h5sclose_f(local_dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
      call h5pclose_f(plist_id, ierr)

   end subroutine hdf5_parallel_write_flt_r2

   !===============================================================================================================================
   subroutine hdf5_parallel_write_int_r1(fid, dset_name, displ, global_dims, idata)

      implicit none
      integer(HID_T), intent(in):: fid
      character(len=*), intent(in):: dset_name
      integer(HSIZE_T), intent(in):: global_dims ! shape of entire dataset being written
      integer(HSSIZE_T), intent(in):: displ ! displacement of process
      integer, intent(in):: idata(:) ! data local to MPI Process to write
      integer(HSIZE_T):: local_dims(1), global_dims_copy(1)
      integer(HSSIZE_T):: displ_copy(1) ! this exists because the hdf5 subroutine only accepts arrays, not constants
      integer, parameter:: rank = 1

      local_dims = size(idata)
      global_dims_copy(1) = global_dims
      displ_copy(1) = displ

      ! Creating dataspace for entire dataset
      call h5screate_simple_f(rank, global_dims_copy, dspace_id, ierr)

      ! Creating dataset for entire dataset
      call h5dcreate_f(fid, dset_name, H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)
      call h5sclose_f(dspace_id, ierr)

      ! Creating dataspace in dataset for each process to write to
      call h5screate_simple_f(rank, local_dims, local_dspace_id, ierr)
      call h5dget_space_f(dset_id, dspace_id, ierr)
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, displ_copy, local_dims, ierr)
      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
      ! Write the dataset collectively.
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, idata, global_dims_copy, ierr, file_space_id=dspace_id, &
                      mem_space_id=local_dspace_id, xfer_prp=plist_id)

      ! closing all dataspaces,datasets,property lists.
      call h5sclose_f(dspace_id, ierr)
      call h5sclose_f(local_dspace_id, ierr)
      call h5dclose_f(dset_id, ierr)
      call h5pclose_f(plist_id, ierr)

   end subroutine hdf5_parallel_write_int_r1
#endif
end module hdf5_parallel_io_helper_m
