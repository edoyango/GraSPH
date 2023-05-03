module output_m

   use datatypes, only: particles
   use hdf5
#ifdef PARALLEL
   use mpi_f08
#endif
   use hdf5_parallel_io_helper_m, only: hdf5_attribute_write, hdf5_parallel_write, hdf5_fileopen_write
   use param, only: f, dim, output_directory, halotype

   private
   public:: output

contains

   !==============================================================================================================================
   subroutine output(itimestep, save_step, my_rank, num_ranks, parts)
      ! Subroutine to write data to disk. Called intermittently, determined by save_step supplied at run-time

      implicit none
      integer, intent(in):: my_rank, num_ranks, itimestep, save_step
      type(particles), intent(in):: parts
      integer:: ntotal_glob(num_ranks), nhalo_glob(num_ranks), nvirt_glob(num_ranks), ierr, n, i, nt, nv, nh
#ifdef PARALLEL
      type(MPI_Request):: request(3)
      type(MPI_Status):: status(3)
#endif
      integer(HID_T):: fid, real_group_id, virt_group_id, halo_group_id, attr_id
      integer(HSIZE_T):: global_dims(2)
      integer(HSSIZE_T):: displ(2)
      character(len=220):: filepath
      character(len=4):: number
      type(particles), allocatable:: real_part_tmp(:), virt_part_tmp(:), halo_part_tmp(:)
      logical, allocatable:: mask(:)

#ifdef PARALLEL
      ! Exchanging how many particles each process will output data for
      call MPI_Iallgather(parts%ntotal_loc, 1, MPI_INTEGER, ntotal_glob, 1, MPI_INTEGER, MPI_COMM_WORLD, request(1), ierr)
      call MPI_Iallgather(parts%nhalo_loc, 1, MPI_INTEGER, nhalo_glob, 1, MPI_INTEGER, MPI_COMM_WORLD, request(2), ierr)
      call MPI_Iallgather(parts%nvirt_loc, 1, MPI_INTEGER, nvirt_glob, 1, MPI_INTEGER, MPI_COMM_WORLD, request(3), ierr)
#else
      ntotal_glob(1) = ntotal_loc
      nvirt_glob(1) = nvirt_loc
      nhalo_glob(1) = 0
#endif

      !Format number
      n = itimestep/save_step
      write (number, '(I4.4)') n

      ! initializing hdf5
      call h5open_f(ierr)

      ! creating hdf5 output file (hdf5_parallel_fileopen is a custom hdf5 wrapper subroutine)
      filepath = trim(output_directory)//"/sph_out"//number//".h5"
      call hdf5_fileopen_write(filepath, fid)

      ! creating groups in hdf5 output for real, halo, virtual, and ghost particles
      call h5gcreate_f(fid, "real", real_group_id, ierr)
      call h5gcreate_f(fid, "virt", virt_group_id, ierr)
      call h5gcreate_f(fid, "halo", halo_group_id, ierr)

      allocate(mask(parts%nsum()))
      
      ! Writing data for real particles ------------------------------------------------------------------------------------------
#ifdef PARALLEL
      ! defining array shapes and displacments
      call MPI_Wait(request(1), status(1), ierr) ! Waiting for non-blocking exchange to complete
#endif

      global_dims(1) = dim
      global_dims(2) = parts%ntotal
      displ(1) = 0 ! hdf5 works with 0-indexing
      if (my_rank==0) then
         displ(2) = 0
      else
         displ(2) = SUM(ntotal_glob(1:my_rank))
      end if

      if (sum(ntotal_glob)/=0) then
         mask = parts%itype(1:parts%nsum()) > 0 .and. parts%itype(1:parts%nsum()) < halotype
         call write_particle_data(my_rank, real_group_id, global_dims, displ, parts, mask, parts%ntotal_loc)
      end if

      call hdf5_attribute_write(real_group_id, parts%ntotal)
      
      call h5gclose_f(real_group_id, ierr)

      ! Writing data for halo particles ------------------------------------------------------------------------------------------
#ifdef PARALLEL
      ! defining array shapes and displacments
      call MPI_WAIT(request(2), status(2), ierr) ! Waiting for non-blocking exchange to complete
#endif

      global_dims(1) = dim
      global_dims(2) = sum(nhalo_glob)
      displ(1) = 0 ! hdf5 works with 0-indexing
      if (my_rank==0) then
         displ(2) = 0
      else
         displ(2) = SUM(nhalo_glob(1:my_rank))
      end if

      if (sum(nhalo_glob)/=0) then
         mask = parts%itype(1:parts%nsum()) < -halotype .or. parts%itype(1:parts%nsum()) > halotype
         call write_particle_data(my_rank, halo_group_id, global_dims, displ, parts, mask, parts%nhalo_loc)
      end if

      call hdf5_attribute_write(halo_group_id, sum(nhalo_glob))

      call h5gclose_f(halo_group_id, ierr)

      ! Writing data for virtual particles ---------------------------------------------------------------------------------------
#ifdef PARALLEL
      ! defining array shapes and displacments
      call MPI_WAIT(request(3), status(3), ierr) ! Waiting for non-blocking exchange to complete
#endif

      global_dims(1) = dim
      global_dims(2) = sum(nvirt_glob)
      displ(1) = 0 ! hdf5 works with 0-indexing
      if (my_rank==0) then
         displ(2) = 0
      else
         displ(2) = SUM(nvirt_glob(1:my_rank))
      end if

      if (sum(nvirt_glob)/=0) then
         mask = parts%itype(1:parts%nsum()) < 0 .and. parts%itype(1:parts%nsum()) > -halotype
         call write_particle_data(my_rank, virt_group_id, global_dims, displ, parts, mask, parts%nvirt_loc)
      end if

      call hdf5_attribute_write(virt_group_id, sum(nvirt_glob))

      call h5gclose_f(virt_group_id, ierr)

      ! Closing output file
      call h5fclose_f(fid, ierr)

      ! closing hdf5
      call h5close_f(ierr)

   end subroutine output

   !==============================================================================================================================
   subroutine write_particle_data(my_rank, gid_in, gdims, ldispl, parts, mask, nelem)

      implicit none
      integer(HID_T), intent(in):: gid_in
      integer, intent(in):: my_rank, nelem
      integer(HSIZE_T), intent(in):: gdims(2)
      integer(HSSIZE_T), intent(in):: ldispl(2)
      type(particles), intent(in):: parts
      logical, intent(in):: mask(:)
      integer:: d, n, i
      integer, allocatable:: int_tmp(:, :)
      real(f), allocatable:: dbl_tmp(:, :)

      allocate (int_tmp(nelem, 1), dbl_tmp(nelem, 1))

      ! writing 1D arrays
      int_tmp(:, 1) = pack(parts%indglob(1:parts%nsum()), mask)
      call hdf5_parallel_write(gid_in, 'ind', ldispl(2), gdims(2), int_tmp(:, 1))
      int_tmp(:, 1) = my_rank
      call hdf5_parallel_write(gid_in, 'procid', ldispl(2), gdims(2), int_tmp(:, 1))
      int_tmp(:, 1) = pack(parts%itype(1:parts%nsum()), mask)
      call hdf5_parallel_write(gid_in, 'type', ldispl(2), gdims(2), int_tmp(:, 1))
      dbl_tmp(:, 1) = pack(parts%rho(1:parts%nsum()), mask)
      call hdf5_parallel_write(gid_in, 'rho', ldispl(2), gdims(2), dbl_tmp(:, 1))
      dbl_tmp(:, 1) = pack(parts%p(1:parts%nsum()), mask)
      call hdf5_parallel_write(gid_in, 'p', ldispl(2), gdims(2), dbl_tmp(:, 1))

      deallocate (int_tmp, dbl_tmp)

      ! Converting position data to contiguous array
      allocate (dbl_tmp(dim, nelem))

      n = 0
      do i = 1, parts%nsum()
         if (mask(i)) then
            n = n + 1
            dbl_tmp(:, n) = parts%x(:, i)
         end if
      end do

      ! writing position
      call hdf5_parallel_write(gid_in, 'x', ldispl, gdims, dbl_tmp)

      ! Converting velocity data to contiguous array
      n = 0
      do i = 1, parts%nsum()
         if (mask(i)) then
            n = n + 1
            dbl_tmp(:, n) = parts%vx(:, i)
         end if
      end do

      ! writing velocity
      call hdf5_parallel_write(gid_in, 'v', ldispl, gdims, dbl_tmp)

      ! cleanup
      deallocate (dbl_tmp)

   end subroutine write_particle_data

end module output_m
