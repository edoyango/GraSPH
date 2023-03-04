module input_m

   use datatypes, only: particles, interactions
   use hdf5
   use param, only: dim, irho, dxo, f, hsml, mass, halotype, input_file
   use hdf5_parallel_io_helper_m, only: hdf5_attribute_read, hdf5_fileopen_read, hdf5_parallel_read
#ifdef PARALLEL
   use mpi
#endif

   public:: read_input_and_allocate, update_virt_part

contains

   !====================================================================================================================
   subroutine allocatePersistentArrays(ntotal, nvirt, parts, pairs, nexti, maxnloc, maxinter)

      implicit none
      integer, intent(in):: ntotal, nvirt
      integer, intent(out), allocatable:: nexti(:)
      type(particles), intent(inout), allocatable, codimension[:]:: parts(:)
      type(interactions), intent(out), allocatable:: pairs(:)
      integer, intent(out):: maxnloc, maxinter

      maxnloc = ntotal + nvirt
      maxinter = 262*maxnloc

      allocate (parts(maxnloc) [*], pairs(maxinter), nexti(maxnloc + 1))

   end subroutine allocatePersistentArrays

   !====================================================================================================================
   subroutine read_input_and_allocate(thisImage, numImages, ntotal, ntotal_loc, nvirt, nvirt_loc, parts, pairs, nexti, &
      maxnloc, maxinter)
      ! Generates initial physical particle configuration.
      ! 2 cases: return only number of particles retrieved, or generating the particles

      implicit none
      integer, intent(in):: thisImage, numImages
      integer, intent(out):: ntotal, ntotal_loc, nvirt, nvirt_loc
      type(particles), allocatable, codimension[:], intent(inout):: parts(:)
      type(interactions), allocatable, intent(out):: pairs(:)
      integer, allocatable, intent(out):: nexti(:)
      integer, intent(out):: maxnloc, maxinter
      integer:: i, j, k, n, ierr, otherImage, ntotal_loc_rounded_avg, nvirt_loc_rounded_avg, nbuff(1)
      logical:: test_init
      integer(HID_T):: fid, gid_r, gid_v, plist_id
      integer(HSIZE_T):: global_dims(2)
      integer(HSSIZE_T):: displ(2)
#ifdef PARALLEL
      ! initializing MPI (if needed)
      call MPI_INITIALIZED(test_init, ierr)
      if (.not. test_init) call MPI_INIT(ierr)
#endif
      ! initializing hdf5
      call h5open_f(ierr)

      ! using helper subroutine to open file
      call hdf5_fileopen_read(input_file, fid)

      ! opening real and virtual particle groups in hdf5 file
      call h5gopen_f(fid, 'real', gid_r, ierr)
      call h5gopen_f(fid, 'virt', gid_v, ierr)

      ! reading the 'n' attribute in each file to get ntotal and nvirt
      call hdf5_attribute_read(gid_r, ntotal)
      call hdf5_attribute_read(gid_v, nvirt)

      ! allocate arrays before reading particle data
      call allocatePersistentArrays(ntotal, nvirt, parts, pairs, nexti, maxnloc, maxinter)

      ! how many particles to read per process. For 1 image, ntotal_loc = ntotal and nvirt_loc = nvirt
      ntotal_loc_rounded_avg = ceiling(dble(ntotal)/numImages)
      nvirt_loc_rounded_avg = ceiling(dble(nvirt)/numImages)
      if (thisImage .eq. numImages) then
         ntotal_loc = ntotal - (numImages - 1)*ntotal_loc_rounded_avg
         nvirt_loc = nvirt - (numImages - 1)*nvirt_loc_rounded_avg
      else
         ntotal_loc = ntotal_loc_rounded_avg
         nvirt_loc = nvirt_loc_rounded_avg
      end if

      ! stopping program if array bounds are exceeded
      !if ((procid .eq. 0) .and. (n_loc .gt. maxnloc)) call error_msg(1, 1)

      ! reading real particle data from hdf5 file
      global_dims(1) = dim
      global_dims(2) = ntotal
      displ(1) = 0
      displ(2) = (thisImage-1)*ntotal_loc_rounded_avg
      if (ntotal /= 0) call read_particle_data_parallel(thisImage, gid_r, global_dims, displ, parts(1:ntotal_loc))

      ! reading virtual particle data from hdf5 file
      global_dims(1) = dim
      global_dims(2) = nvirt
      displ(1) = 0
      displ(2) = (thisImage-1)*nvirt_loc_rounded_avg
      if (nvirt /= 0) call read_particle_data_parallel(thisImage, gid_v, global_dims, displ, &
         parts(ntotal_loc+1:ntotal_loc+nvirt_loc))

      do i = 1, ntotal_loc+nvirt_loc
         parts(i)%indloc = i
      end do

      ! closing groups and hdf5 file
      call h5gclose_f(gid_r, ierr)
      call h5gclose_f(gid_v, ierr)
      call h5fclose_f(fid, ierr)

      ! closing hdf5
      call h5close_f(ierr)

#ifdef PARALLEL
      if (.not. test_init) call MPI_FINALIZE(ierr)
#endif

   end subroutine read_input_and_allocate

   !==============================================================================================================================
   subroutine update_virt_part(ntotal_loc, nhalo_loc, nvirt_loc, parts, niac, pairs, nexti, vw)

      implicit none
      integer, intent(in):: ntotal_loc, nhalo_loc, nvirt_loc, niac, nexti(:)
      type(interactions), intent(in):: pairs(:)
      type(particles), intent(inout):: parts(:)
      real(f), intent(inout):: vw(:)
      integer:: i, j, k, tmptype
      real(f):: tmp

      vw(:) = 0._f

      do i = 1, ntotal_loc + nvirt_loc + nhalo_loc
         if (parts(i)%itype<0) then
            parts(i)%rho = 0._f
            parts(i)%vx(:) = 0._f
         end if
      end do

      do i = 1, ntotal_loc + nhalo_loc + nvirt_loc
         do k = nexti(i), nexti(i + 1) - 1
            j = pairs(k)%j
            if ( parts(i)%itype < 0 .and. parts(j)%itype > 0) then
               if (parts(i)%itype < -halotype) then
                  tmptype = parts(i)%itype + halotype
               else
                  tmptype = parts(i)%itype
               end if
               tmp = mass*pairs(k)%w/parts(j)%rho
               vw(i) = vw(i) + tmp
               parts(i)%rho = parts(i)%rho + mass*pairs(k)%w
               select case (tmptype)
               case default
                  parts(i)%vx(:) = parts(i)%vx(:) - parts(j)%vx(:)*tmp
               case (-2)
                  parts(i)%vx(1) = parts(i)%vx(1) + parts(j)%vx(1)*tmp
                  parts(i)%vx(2) = parts(i)%vx(2) + parts(j)%vx(2)*tmp
                  parts(i)%vx(3) = parts(i)%vx(3) - parts(j)%vx(3)*tmp
               case (-3)
                  parts(i)%vx(1) = parts(i)%vx(1) - parts(j)%vx(1)*tmp
                  parts(i)%vx(2) = parts(i)%vx(2) + parts(j)%vx(2)*tmp
                  parts(i)%vx(3) = parts(i)%vx(3) + parts(j)%vx(3)*tmp
               case (-4)
                  parts(i)%vx(1) = parts(i)%vx(1) + parts(j)%vx(1)*tmp
                  parts(i)%vx(2) = parts(i)%vx(2) - parts(j)%vx(2)*tmp
                  parts(i)%vx(3) = parts(i)%vx(3) + parts(j)%vx(3)*tmp
               case (-5)
                  parts(i)%vx(1) = parts(i)%vx(1) - parts(j)%vx(1)*tmp
                  parts(i)%vx(2) = parts(i)%vx(2) - parts(j)%vx(2)*tmp
                  parts(i)%vx(3) = parts(i)%vx(3) + parts(j)%vx(3)*tmp
               end select
            else if (parts(j)%itype < 0 .and. parts(i)%itype > 0) then
               if (parts(j)%itype < -halotype) then
                  tmptype = parts(j)%itype + halotype
               else
                  tmptype = parts(j)%itype
               end if
               tmp = mass*pairs(k)%w/parts(i)%rho
               vw(j) = vw(j) + tmp
               parts(j)%rho = parts(j)%rho + mass*pairs(k)%w
               select case (tmptype)
               case default
                  parts(j)%vx(:) = parts(j)%vx(:) - parts(i)%vx(:)*tmp
               case (-2)
                  parts(j)%vx(1) = parts(j)%vx(1) + parts(i)%vx(1)*tmp
                  parts(j)%vx(2) = parts(j)%vx(2) + parts(i)%vx(2)*tmp
                  parts(j)%vx(3) = parts(j)%vx(3) - parts(i)%vx(3)*tmp
               case (-3)
                  parts(j)%vx(1) = parts(j)%vx(1) - parts(i)%vx(1)*tmp
                  parts(j)%vx(2) = parts(j)%vx(2) + parts(i)%vx(2)*tmp
                  parts(j)%vx(3) = parts(j)%vx(3) + parts(i)%vx(3)*tmp
               case (-4)
                  parts(j)%vx(1) = parts(j)%vx(1) + parts(i)%vx(1)*tmp
                  parts(j)%vx(2) = parts(j)%vx(2) - parts(i)%vx(2)*tmp
                  parts(j)%vx(3) = parts(j)%vx(3) + parts(i)%vx(3)*tmp
               case (-5)
                  parts(j)%vx(1) = parts(j)%vx(1) - parts(i)%vx(1)*tmp
                  parts(j)%vx(2) = parts(j)%vx(2) - parts(i)%vx(2)*tmp
                  parts(j)%vx(3) = parts(j)%vx(3) + parts(i)%vx(3)*tmp
               end select
            end if
         end do
      end do

      do i = 1, ntotal_loc+nvirt_loc+nhalo_loc
         if (parts(i)%itype<0) then
            if (vw(i) > 0._f) then
               parts(i)%rho = parts(i)%rho/vw(i)
               parts(i)%vx(:) = parts(i)%vx(:)/vw(i)
            else
               parts(i)%rho = irho
               parts(i)%vx(:) = 0._f
            end if
         end if
      end do

   end subroutine update_virt_part

   !====================================================================================================================
   subroutine read_particle_data_parallel(thisImage, gid_in, gdims, ldispl, parts)

      implicit none
      integer(HID_T), intent(in):: gid_in
      integer, intent(in):: thisImage
      integer(HSIZE_T), intent(in):: gdims(2)
      integer(HSSIZE_T), intent(in):: ldispl(2)
      type(particles), intent(inout):: parts(:)
      integer:: nelem, i
      integer, allocatable:: int_tmp(:, :)
      real(f), allocatable:: dbl_tmp(:, :)

      nelem = size(parts, 1)
      allocate(int_tmp(nelem, 1), dbl_tmp(nelem, 1))

      ! read 1d arrays
      call hdf5_parallel_read(gid_in, 'ind', ldispl(2), gdims(2), int_tmp(:, 1))
      parts(:)%indglob = int_tmp(:, 1)
      
      call hdf5_parallel_read(gid_in, 'type', ldispl(2), gdims(2), int_tmp(:, 1))
      parts(:)%itype = int_tmp(:, 1)

      call hdf5_parallel_read(gid_in, 'rho', ldispl(2), gdims(2), dbl_tmp(:, 1))
      parts(:)%rho = dbl_tmp(:, 1)

      call hdf5_parallel_read(gid_in, 'p', ldispl(2), gdims(2), dbl_tmp(:, 1))
      parts(:)%p = dbl_tmp(:, 1)

      deallocate(int_tmp, dbl_tmp)

      allocate( dbl_tmp(dim,nelem))
      
      call hdf5_parallel_read(gid_in, 'x', ldispl, gdims, dbl_tmp)
      do i = 1, nelem
         parts(i)%x(:) = dbl_tmp(:, i)
      end do

      call hdf5_parallel_read(gid_in, 'v', ldispl, gdims, dbl_tmp)
      do i = 1, nelem
         parts(i)%vx(:) = dbl_tmp(:, i)
      end do

   end subroutine read_particle_data_parallel

end module input_m
