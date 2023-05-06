module input_m

   use datatypes, only: particles, interactions
   use hdf5
   use param, only: dim, irho, dxo, f, hsml, mass, halotype, input_file
   use hdf5_parallel_io_helper_m, only: hdf5_attribute_read, hdf5_fileopen_read, hdf5_parallel_read
#ifdef PARALLEL
   use mpi
#endif

   public:: read_input_and_allocate!, update_virt_part

contains

   !====================================================================================================================
   subroutine allocatePersistentArrays(parts, pairs, nexti, maxinter)

      implicit none
      integer, intent(out), allocatable:: nexti(:)
      type(particles), intent(inout):: parts
      type(interactions), intent(out), allocatable:: pairs(:)
      integer, intent(out):: maxinter

      parts%maxn = parts%ntotal + parts%nvirt
      maxinter = 262*parts%maxn

      call parts%allocate_particles()

      allocate (pairs(maxinter), nexti(parts%maxn + 1))

   end subroutine allocatePersistentArrays

   !====================================================================================================================
   subroutine read_input_and_allocate(my_rank, num_ranks, parts, pairs, nexti, maxinter)
      ! Generates initial physical particle configuration.
      ! 2 cases: return only number of particles retrieved, or generating the particles

      implicit none
      integer, intent(in):: my_rank, num_ranks
      type(particles), intent(inout):: parts
      type(interactions), allocatable, intent(out):: pairs(:)
      integer, allocatable, intent(out):: nexti(:)
      integer, intent(out):: maxinter
      integer:: i, j, k, n, ierr, otherImage, ntotal_loc_rounded_avg, nvirt_loc_rounded_avg, nbuff(1)
      logical:: test_init
      integer(HID_T):: fid, gid_r, gid_v, plist_id
      integer(HSIZE_T):: global_dims(2)
      integer(HSSIZE_T):: displ(2)

      ! initializing hdf5
      call h5open_f(ierr)

      ! using helper subroutine to open file
      call hdf5_fileopen_read(input_file, fid)

      ! opening real and virtual particle groups in hdf5 file
      call h5gopen_f(fid, 'real', gid_r, ierr)
      call h5gopen_f(fid, 'virt', gid_v, ierr)

      ! reading the 'n' attribute in each file to get ntotal and nvirt
      call hdf5_attribute_read(gid_r, parts%ntotal)
      call hdf5_attribute_read(gid_v, parts%nvirt)

      ! allocate arrays before reading particle data
      call allocatePersistentArrays(parts, pairs, nexti, maxinter)

      ! how many particles to read per process. For 1 image, ntotal_loc = ntotal and nvirt_loc = nvirt
      ntotal_loc_rounded_avg = ceiling(real(parts%ntotal)/num_ranks)
      nvirt_loc_rounded_avg = ceiling(real(parts%nvirt)/num_ranks)
      if (my_rank == num_ranks-1) then
         parts%ntotal_loc = parts%ntotal - (num_ranks - 1)*ntotal_loc_rounded_avg
         parts%nvirt_loc = parts%nvirt - (num_ranks - 1)*nvirt_loc_rounded_avg
      else
         parts%ntotal_loc = ntotal_loc_rounded_avg
         parts%nvirt_loc = nvirt_loc_rounded_avg
      end if

      ! stopping program if array bounds are exceeded
      !if ((procid .eq. 0) .and. (n_loc .gt. maxnloc)) call error_msg(1, 1)

      ! reading real particle data from hdf5 file
      global_dims(1) = dim
      global_dims(2) = parts%ntotal
      displ(1) = 0
      displ(2) = my_rank*ntotal_loc_rounded_avg
      if (parts%ntotal /= 0) call read_particle_data_parallel(my_rank, gid_r, global_dims, displ, parts, 1, &
         parts%ntotal_loc)

      ! reading virtual particle data from hdf5 file
      global_dims(1) = dim
      global_dims(2) = parts%nvirt
      displ(1) = 0
      displ(2) = my_rank*nvirt_loc_rounded_avg
      if (parts%nvirt /= 0) call read_particle_data_parallel(my_rank, gid_v, global_dims, displ, parts, &
         parts%ntotal_loc+1, parts%ntotal_loc+parts%nvirt_loc)

      do i = 1, parts%ntotal_loc+parts%nvirt_loc
         parts%indloc(i) = i
      end do

      ! closing groups and hdf5 file
      call h5gclose_f(gid_r, ierr)
      call h5gclose_f(gid_v, ierr)
      call h5fclose_f(fid, ierr)

      ! closing hdf5
      call h5close_f(ierr)

   end subroutine read_input_and_allocate

   !==============================================================================================================================
   subroutine update_virt_part(parts, niac, pairs, nexti, vw)

      use kernel_m, only: kernel_w

      implicit none
      integer, intent(in):: niac, nexti(:)
      type(interactions), intent(in):: pairs(:)
      type(particles), intent(inout):: parts
      real(f), intent(inout):: vw(:)
      integer:: i, j, k, tmptype
      real(f):: tmp, tw

      do i = 1, parts%nsum()
         if (parts%itype(i)<0) then
            parts%rho(i) = 0._f
            parts%vx(:,i) = 0._f
            vw(i) = 0._f
         end if
      end do

      do i = 1, parts%nsum()
         do k = nexti(i), nexti(i + 1) - 1
            j = pairs(k)%j
            if ( parts%itype(i) < 0 .and. parts%itype(j) > 0) then
               if (parts%itype(i) < -halotype) then
                  tmptype = parts%itype(i) + halotype
               else
                  tmptype = parts%itype(i)
               end if
               tw = kernel_w(sqrt(sum(pairs(k)%dx*pairs(k)%dx)), hsml)
               tmp = mass*tw/parts%rho(j)
               vw(i) = vw(i) + tmp
               parts%rho(i) = parts%rho(i) + mass*tw
               select case (tmptype)
               case default
                  parts%vx(:, i) = parts%vx(:, i) - parts%vx(:, j)*tmp
               case (-2)
                  parts%vx(1, i) = parts%vx(1, i) + parts%vx(1, j)*tmp
                  parts%vx(2, i) = parts%vx(2, i) + parts%vx(2, j)*tmp
                  parts%vx(3, i) = parts%vx(3, i) - parts%vx(3, j)*tmp
               case (-3)
                  parts%vx(1, i) = parts%vx(1, i) - parts%vx(1, j)*tmp
                  parts%vx(2, i) = parts%vx(2, i) + parts%vx(2, j)*tmp
                  parts%vx(3, i) = parts%vx(3, i) + parts%vx(3, j)*tmp
               case (-4)
                  parts%vx(1, i) = parts%vx(1, i) + parts%vx(1, j)*tmp
                  parts%vx(2, i) = parts%vx(2, i) - parts%vx(2, j)*tmp
                  parts%vx(3, i) = parts%vx(3, i) + parts%vx(3, j)*tmp
               case (-5)
                  parts%vx(1, i) = parts%vx(1, i) - parts%vx(1, j)*tmp
                  parts%vx(2, i) = parts%vx(2, i) - parts%vx(2, j)*tmp
                  parts%vx(3, i) = parts%vx(3, i) + parts%vx(3, j)*tmp
               end select
            else if (parts%itype(j) < 0 .and. parts%itype(i) > 0) then
               if (parts%itype(j) < -halotype) then
                  tmptype = parts%itype(j) + halotype
               else
                  tmptype = parts%itype(j)
               end if
               tw = kernel_w(sqrt(sum(pairs(k)%dx*pairs(k)%dx)), hsml)
               tmp = mass*tw/parts%rho(i)
               vw(j) = vw(j) + tmp
               parts%rho(j) = parts%rho(j) + mass*tw
               select case (tmptype)
               case default
                  parts%vx(:, j) = parts%vx(:, j) - parts%vx(:, i)*tmp
               case (-2)
                  parts%vx(1, j) = parts%vx(1, j) + parts%vx(1, i)*tmp
                  parts%vx(2, j) = parts%vx(2, j) + parts%vx(2, i)*tmp
                  parts%vx(3, j) = parts%vx(3, j) - parts%vx(3, i)*tmp
               case (-3)
                  parts%vx(1, j) = parts%vx(1, j) - parts%vx(1, i)*tmp
                  parts%vx(2, j) = parts%vx(2, j) + parts%vx(2, i)*tmp
                  parts%vx(3, j) = parts%vx(3, j) + parts%vx(3, i)*tmp
               case (-4)
                  parts%vx(1, j) = parts%vx(1, j) + parts%vx(1, i)*tmp
                  parts%vx(2, j) = parts%vx(2, j) - parts%vx(2, i)*tmp
                  parts%vx(3, j) = parts%vx(3, j) + parts%vx(3, i)*tmp
               case (-5)
                  parts%vx(1, j) = parts%vx(1, j) - parts%vx(1, i)*tmp
                  parts%vx(2, j) = parts%vx(2, j) - parts%vx(2, i)*tmp
                  parts%vx(3, j) = parts%vx(3, j) + parts%vx(3, i)*tmp
               end select
            end if
         end do
      end do

      do i = 1, parts%nsum()
         if (parts%itype(i)<0) then
            if (vw(i) > 0._f) then
               parts%rho(i) = parts%rho(i)/vw(i)
               parts%vx(:, i) = parts%vx(:, i)/vw(i)
            else
               parts%rho(i) = irho
               parts%vx(:, i) = 0._f
            end if
         end if
      end do

   end subroutine update_virt_part

   !====================================================================================================================
   subroutine read_particle_data_parallel(my_rank, gid_in, gdims, ldispl, parts, istart, iend)

      implicit none
      integer(HID_T), intent(in):: gid_in
      integer, intent(in):: my_rank, istart, iend
      integer(HSIZE_T), intent(in):: gdims(2)
      integer(HSSIZE_T), intent(in):: ldispl(2)
      type(particles), intent(inout):: parts
      integer:: i

      ! read 1d arrays
      call hdf5_parallel_read(gid_in, 'ind', ldispl(2), gdims(2), parts%indglob(istart:iend))
      call hdf5_parallel_read(gid_in, 'type', ldispl(2), gdims(2), parts%itype(istart:iend))
      call hdf5_parallel_read(gid_in, 'rho', ldispl(2), gdims(2), parts%rho(istart:iend))
      call hdf5_parallel_read(gid_in, 'p', ldispl(2), gdims(2), parts%p(istart:iend))
      
      ! read 2d arrays
      call hdf5_parallel_read(gid_in, 'x', ldispl, gdims, parts%x(:, istart:iend))
      call hdf5_parallel_read(gid_in, 'v', ldispl, gdims, parts%vx(:, istart:iend))

   end subroutine read_particle_data_parallel

end module input_m
