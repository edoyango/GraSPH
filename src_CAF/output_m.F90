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
   subroutine output(itimestep, save_step, thisImage, numImages, ntotal_loc, nhalo_loc, nvirt_loc, parts, ntotal)
      ! Subroutine to write data to disk. Called intermittently, determined by save_step supplied at run-time

      implicit none
      integer, intent(in):: thisImage, numImages, itimestep, save_step, ntotal_loc, nhalo_loc, nvirt_loc, ntotal
      type(particles), intent(in):: parts(:)
      integer:: ntotal_glob(numImages), nhalo_glob(numImages), nvirt_glob(numImages), ierr, n, i, nt, nv, nh
#ifdef PARALLEL
      type(MPI_Request):: request(3)
      type(MPI_Status):: status(3)
#endif
      integer(HID_T):: fid, real_group_id, virt_group_id, halo_group_id, attr_id
      integer(HSIZE_T):: global_dims(2)
      integer(HSSIZE_T):: displ(2)
      character(len=220):: filepath
      character(len=4):: number
      logical:: test_init
      type(particles), allocatable:: real_part_tmp(:), virt_part_tmp(:), halo_part_tmp(:)

#ifdef PARALLEL
      call MPI_INITIALIZED(test_init, ierr)
      if (.not. test_init) call MPI_INIT(ierr)

      ! Exchanging how many particles each process will output data for
      call MPI_IALLGATHER(ntotal_loc, 1, MPI_INTEGER, ntotal_glob, 1, MPI_INTEGER, MPI_COMM_WORLD, request(1), ierr)
      call MPI_IALLGATHER(nhalo_loc, 1, MPI_INTEGER, nhalo_glob, 1, MPI_INTEGER, MPI_COMM_WORLD, request(2), ierr)
      call MPI_IALLGATHER(nvirt_loc, 1, MPI_INTEGER, nvirt_glob, 1, MPI_INTEGER, MPI_COMM_WORLD, request(3), ierr)
#else
      ntotal_glob(1) = ntotal_loc
      nvirt_glob(1) = nvirt_loc
      nhalo_glob(1) = 0
#endif

      ! putting real, virtual, halo particles into temporary arrays
      ! may look redundant, but is needed because particles are not guaranteed to be ordered by type
      allocate(real_part_tmp(ntotal_loc), virt_part_tmp(nvirt_loc), halo_part_tmp(nhalo_loc))

      nt = 0
      nv = 0
      nh = 0
      do i = 1, ntotal_loc+nvirt_loc+nhalo_loc
         if (parts(i)%itype>0 .and. parts(i)%itype<halotype) then
            nt = nt + 1
            real_part_tmp(nt) = parts(i)
         else if (parts(i)%itype<0 .and. parts(i)%itype>-halotype) then
            nv = nv + 1
            virt_part_tmp(nv) = parts(i)
         else if (parts(i)%itype>halotype .or. parts(i)%itype<-halotype) then
            nh = nh + 1
            halo_part_tmp(nh) = parts(i)
         end if
      end do

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
      
      ! Writing data for real particles ------------------------------------------------------------------------------------------
#ifdef PARALLEL
      ! defining array shapes and displacments
      call MPI_WAIT(request(1), status(1), ierr) ! Waiting for non-blocking exchange to complete
#endif

      global_dims(1) = dim
      global_dims(2) = ntotal
      displ(1) = 0 ! hdf5 works with 0-indexing
      if (thisimage==1) then
         displ(2) = 0
      else
         displ(2) = SUM(ntotal_glob(1:thisImage - 1))
      end if

      if (sum(ntotal_glob)/=0) then
         call write_particle_data(thisImage, real_group_id, global_dims, displ, real_part_tmp)
      end if

      call hdf5_attribute_write(real_group_id, ntotal)
      
      call h5gclose_f(real_group_id, ierr)

      ! Writing data for halo particles ------------------------------------------------------------------------------------------
#ifdef PARALLEL
      ! defining array shapes and displacments
      call MPI_WAIT(request(2), status(2), ierr) ! Waiting for non-blocking exchange to complete
#endif

      global_dims(1) = dim
      global_dims(2) = sum(nhalo_glob)
      displ(1) = 0 ! hdf5 works with 0-indexing
      if (thisimage==1) then
         displ(2) = 0
      else
         displ(2) = SUM(nhalo_glob(1:thisImage - 1))
      end if

      if (sum(nhalo_glob)/=0) then
         call write_particle_data(thisImage, halo_group_id, global_dims, displ, halo_part_tmp)
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
      if (thisimage==1) then
         displ(2) = 0
      else
         displ(2) = SUM(nvirt_glob(1:thisImage - 1))
      end if

      if (sum(nvirt_glob)/=0) then
         call write_particle_data(thisImage, virt_group_id, global_dims, displ, virt_part_tmp)
      end if

      call hdf5_attribute_write(virt_group_id, sum(nvirt_glob))

      call h5gclose_f(virt_group_id, ierr)

      ! Closing output file
      call h5fclose_f(fid, ierr)

      ! closing hdf5
      call h5close_f(ierr)
#ifdef PARALLEL
      if (.not. test_init) call MPI_FINALIZE(ierr)
#endif

   end subroutine output

   !==============================================================================================================================
   subroutine write_particle_data(thisImage, gid_in, gdims, ldispl, parts)

      implicit none
      integer(HID_T), intent(in):: gid_in
      integer, intent(in):: thisImage
      integer(HSIZE_T), intent(in):: gdims(2)
      integer(HSSIZE_T), intent(in):: ldispl(2)
      type(particles), intent(in):: parts(:)
      integer:: nelem, i
      integer, allocatable:: int_tmp(:, :)
      real(f), allocatable:: dbl_tmp(:, :)

      nelem = size(parts, 1)
      allocate (int_tmp(nelem, 1), dbl_tmp(nelem, 1))

      ! writing 1D arrays
      int_tmp(:, 1) = parts(:)%indglob
      call hdf5_parallel_write(gid_in, 'ind', ldispl(2), gdims(2), int_tmp(:, 1))
      int_tmp(:, 1) = thisImage
      call hdf5_parallel_write(gid_in, 'procid', ldispl(2), gdims(2), int_tmp(:, 1))
      int_tmp(:, 1) = parts(:)%itype
      call hdf5_parallel_write(gid_in, 'type', ldispl(2), gdims(2), int_tmp(:, 1))
      dbl_tmp(:, 1) = parts(:)%rho
      call hdf5_parallel_write(gid_in, 'rho', ldispl(2), gdims(2), dbl_tmp(:, 1))
      dbl_tmp(:, 1) = parts(:)%p
      call hdf5_parallel_write(gid_in, 'p', ldispl(2), gdims(2), dbl_tmp(:, 1))

      deallocate (int_tmp, dbl_tmp)

      ! Converting position data to contiguous array
      allocate (dbl_tmp(dim, nelem))

      do i = 1, nelem
         dbl_tmp(:, i) = parts(i)%x(:)
      end do

      ! writing position
      call hdf5_parallel_write(gid_in, 'x', ldispl, gdims, dbl_tmp)

      ! Converting velocity data to contiguous array
      do i = 1, nelem
         dbl_tmp(:, i) = parts(i)%vx(:)
      end do

      ! writing velocity
      call hdf5_parallel_write(gid_in, 'v', ldispl, gdims, dbl_tmp)

      ! cleanup
      deallocate (dbl_tmp)

   end subroutine write_particle_data

end module output_m
