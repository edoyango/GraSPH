module output_m
   
   use datatypes, only: particles
   use param, only: f, dim, output_directory, output_phys, output_halo, output_virt, output_flt_type

   use hdf5
   use mpi_f08

   use hdf5_parallel_io_helper_m, only: hdf5_parallel_fileopen, hdf5_parallel_write

   private
   public:: output

contains

   !==============================================================================================================================
   subroutine output(itimestep, save_step, procid, numprocs, ntotal_loc, nhalo_loc, nvirt_loc, nghos_loc, parts, ntotal)
      ! Subroutine to write data to disk. Called intermittently, determined by save_step supplied at run-time

      implicit none
      integer, intent(in):: procid, numprocs, itimestep, save_step, ntotal_loc, nhalo_loc, nvirt_loc, nghos_loc, ntotal
      type(particles), intent(in):: parts(:)
      integer:: ntotal_glob(numprocs), nhalo_glob(numprocs), nvirt_glob(numprocs), nghos_glob(numprocs), ierr, n, i
      type(MPI_Request):: request(4)
      type(MPI_Status):: status(4)
      integer(HID_T):: fid, real_group_id, virt_group_id, ghos_group_id
      integer(HSIZE_T):: global_dims(2)
      integer(HSSIZE_T):: displ(2)
      character(len=220):: filepath
      character(len=4):: number
      
      ! Exchanging how many particles each process will output data for
      call MPI_IALLGATHER(ntotal_loc, 1, MPI_INTEGER, ntotal_glob, 1, MPI_INTEGER, MPI_COMM_WORLD, request(1), ierr)
      call MPI_IALLGATHER(nhalo_loc, 1, MPI_INTEGER, nhalo_glob, 1, MPI_INTEGER, MPI_COMM_WORLD, request(2), ierr)
      call MPI_IALLGATHER(nvirt_loc, 1, MPI_INTEGER, nvirt_glob, 1, MPI_INTEGER, MPI_COMM_WORLD, request(3), ierr)
      call MPI_IALLGATHER(nghos_loc, 1, MPI_INTEGER, nghos_glob, 1, MPI_INTEGER, MPI_COMM_WORLD, request(4), ierr)

      !Format number
      n = itimestep/save_step
      write (number, '(I4.4)') n

      ! initializing hdf5
      call h5open_f(ierr)

      ! creating hdf5 output file (hdf5_parallel_fileopen is a custom hdf5 wrapper subroutine)
      filepath = trim(output_directory)//"/sph_out"//number//".h5"
      call hdf5_parallel_fileopen(filepath, fid)

      ! creating groups in hdf5 output for real, halo, virtual, and ghost particles
      call h5gcreate_f(fid, "real", real_group_id, ierr)
      call h5gclose_f(real_group_id, ierr)

      call h5gcreate_f(fid, "virt", virt_group_id, ierr)
      call h5gclose_f(virt_group_id, ierr)

      call h5gcreate_f(fid, "halo", ghos_group_id, ierr)
      call h5gclose_f(ghos_group_id, ierr)

      call h5gcreate_f(fid, "ghos", ghos_group_id, ierr)
      call h5gclose_f(ghos_group_id, ierr)

      ! Writing data for real particles ------------------------------------------------------------------------------------------
      ! defining array shapes and displacments
      call MPI_WAIT(request(1), status(1), ierr) ! Waiting for non-blocking exchange to complete

      global_dims(:) = [dim, ntotal]
      displ(:) = [0, SUM(ntotal_glob(1:procid))] ! hdf5 works with 0-indexing

      call write_particle_data(procid, fid, 'real', global_dims, displ, parts(1:ntotal_loc))

      ! Writing data for halo particles ------------------------------------------------------------------------------------------
      ! defining array shapes and displacments
      call MPI_WAIT(request(2), status(2), ierr) ! Waiting for non-blocking exchange to complete

      global_dims(:) = [dim, sum(nhalo_glob)]
      displ(:) = [0, SUM(nhalo_glob(1:procid))] ! hdf5 works with 0-indexing

      call write_particle_data(procid, fid, 'halo', global_dims, displ, parts(ntotal_loc+1:ntotal_loc+nhalo_loc))

      ! Writing data for virtual particles ---------------------------------------------------------------------------------------
      ! defining array shapes and displacments
      call MPI_WAIT(request(3), status(3), ierr) ! Waiting for non-blocking exchange to complete

      global_dims(:) = [dim, sum(nvirt_glob)]
      displ(:) = [0, SUM(nvirt_glob(1:procid))] ! hdf5 works with 0-indexing

      call write_particle_data(procid, fid, 'virt', global_dims, displ, &
      parts(ntotal_loc + nhalo_loc + 1: ntotal_loc + nhalo_loc + nvirt_loc))

      ! Writing data for ghost particles -----------------------------------------------------------------------------------------
      ! defining array shapes and displacments
      call MPI_WAIT(request(4), status(4), ierr) ! Waiting for non-blocking exchange to complete

      global_dims(:) = [dim, sum(nghos_glob)]
      displ(:) = [0, SUM(nghos_glob(1:procid))] ! hdf5 works with 0-indexing

      call write_particle_data(procid, fid, 'ghos', global_dims, displ, &
      parts(ntotal_loc + nhalo_loc + nvirt_loc +1: ntotal_loc + nhalo_loc + nvirt_loc + nghos_loc))

      ! Closing output file
      call h5fclose_f(fid, ierr)

      ! closing hdf5
      call h5close_f(ierr)

   end subroutine output

   !==============================================================================================================================
   subroutine write_particle_data(procid, fid_in, particletype, gdims, ldispl, parts)

      implicit none
      integer(HID_T), intent(in):: fid_in
      character(*), intent(in):: particletype
      integer, intent(in):: procid
      integer(HSIZE_T), intent(in):: gdims(2)
      integer(HSSIZE_T), intent(in):: ldispl(2)
      type(particles), intent(in):: parts(:)
      integer:: nelem, i
      integer, allocatable:: int_tmp(:, :)
      real(f), allocatable:: dbl_tmp(:, :)

      nelem = size(parts,1)
      allocate (int_tmp(nelem, 1), dbl_tmp(nelem, 1))

      ! writing 1D arrays
      int_tmp(:, 1) = parts(:)%indglob
      call hdf5_parallel_write(fid_in, particletype//'/ind', ldispl(2), gdims(2), int_tmp(:, 1))
      int_tmp(:, 1) = procid
      call hdf5_parallel_write(fid_in, particletype//'/procid', ldispl(2), gdims(2), int_tmp(:, 1))
      int_tmp(:, 1) = parts(:)%itype
      call hdf5_parallel_write(fid_in, particletype//'/type', ldispl(2), gdims(2), int_tmp(:, 1))
      dbl_tmp(:, 1) = parts(:)%rho
      call hdf5_parallel_write(fid_in, particletype//'/rho', ldispl(2), gdims(2), dbl_tmp(:, 1))
      dbl_tmp(:, 1) = parts(:)%p
      call hdf5_parallel_write(fid_in, particletype//'/p', ldispl(2), gdims(2), dbl_tmp(:, 1))

      deallocate (int_tmp, dbl_tmp)

      ! Converting position data to contiguous array
      allocate (dbl_tmp(dim, nelem))

      do i = 1, nelem
         dbl_tmp(:, i) = parts(i)%x(:)
      end do

      ! writing position
      call hdf5_parallel_write(fid_in, particletype//'/x', ldispl, gdims, dbl_tmp)

      ! Converting velocity data to contiguous array
      do i = 1, nelem
         dbl_tmp(:, i) = parts(i)%vx(:)
      end do

      ! writing velocity
      call hdf5_parallel_write(fid_in, particletype//'/v', ldispl, gdims, dbl_tmp)

      ! cleanup
      deallocate (dbl_tmp)

   end subroutine write_particle_data

end module output_m
