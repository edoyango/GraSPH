module input_m

   use datatypes, only: particles, interactions
   use hdf5
   use param, only: dim, irho, dxo, f, hsml, mp, np, op, pp, qp, rp, nlayer, mass, halotype, input_file
   use hdf5_parallel_io_helper_m, only: hdf5_parallel_fileopen_read, hdf5_parallel_read
   use mpi

   private
   real(f), parameter:: vxmin = 0._f, vymin = 0._f, vzmin = 0._f, &
                        vxmax = vxmin + pp*dxo, vymax = vymin + qp*dxo, vzmax = vzmin + rp*dxo
   real(f), parameter:: rxmin = 0._f, rymin = 0._f, rzmin = 0._f, &
                        rxmax = rxmin + mp*dxo, rymax = rymin + np*dxo, rzmax = rzmin + op*dxo

   public:: return_ntotal, return_nvirt, allocatePersistentArrays, generate_real_part, generate_virt_part, &
            update_virt_part

contains

   !====================================================================================================================
   function return_ntotal() result(ntotal)

      implicit none
      integer:: ntotal

      ntotal = mp*np*op

   end function

   !====================================================================================================================
   function return_nvirt() result(nvirt)

      implicit none
      integer:: nvirt

      nvirt = (rp + 2*nlayer)*(2*nlayer + pp)*(2*nlayer + qp) - pp*qp*rp

   end function

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
   subroutine generate_real_part(thisImage, numImages, ntotal, ntotal_loc, parts)
      ! Generates initial physical particle configuration.
      ! 2 cases: return only number of particles retrieved, or generating the particles

      implicit none
      integer, intent(in):: thisImage, numImages, ntotal
      integer, intent(out):: ntotal_loc
      type(particles), intent(out):: parts(:)
      integer:: i, j, k, n, ierr, otherImage
      integer, allocatable, codimension[:]:: ntotal_glob(:)
      logical:: test_init
      integer(HID_T):: fid
      integer(HSIZE_T):: global_dims(2)
      integer(HSSIZE_T):: displ(2)

      ! how many particles to generate per process
      ntotal_loc = ceiling(dble(ntotal)/numImages)
      if (thisImage .eq. numImages) ntotal_loc = ntotal - (numImages - 1)*ntotal_loc

      ! stopping program if array bounds are exceeded
      !if ((procid .eq. 0) .and. (n_loc .gt. maxnloc)) call error_msg(1, 1)

      ! Exchanging how many particles each process will read data for
      allocate(ntotal_glob(numImages)[*])
      do n = 0, numImages-1
         otherImage = mod(thisImage+n, numImages)
         if (otherImage==0) otherImage = numImages
         ntotal_glob(thisImage) [otherImage] = ntotal_loc
      end do

      call MPI_INITIALIZED(test_init, ierr)
      if (.not. test_init) call MPI_INIT(ierr)

      ! initializing hdf5
      call h5open_f(ierr)

      call hdf5_parallel_fileopen_read(input_file, fid)

      sync all

      global_dims(1) = dim; global_dims(2) = ntotal
      displ(1) = 0; displ(2) = SUM(ntotal_glob(1:thisImage-1))
      if (ntotal /= 0) call read_particle_data_parallel(thisImage, fid, 'real', global_dims, displ, parts(1:ntotal_loc))

      do i = 1, ntotal_loc
         parts(i)%indloc = i
      end do

      ! Closing output file
      call h5fclose_f(fid, ierr)

      ! closing hdf5
      call h5close_f(ierr)

      if (.not. test_init) call MPI_FINALIZE(ierr)

   end subroutine generate_real_part

   !====================================================================================================================
   subroutine generate_virt_part(thisImage, numImages, ntotal, ntotal_loc, nvirt, nvirt_loc, parts)

      ! Generates the virtual particle configuration. Can change over time or remain static
      ! 2 cases: return only number of particles retrieved, or generating the particles

      implicit none
      integer, intent(in):: thisImage, numImages, ntotal_loc, ntotal, nvirt
      type(particles), intent(inout):: parts(:)
      integer, intent(out):: nvirt_loc
      integer:: i, j, k, n, n_loc, n_loc_i, n_start, n_done, ierr, otherImage
      integer, allocatable, codimension[:]:: nvirt_glob(:)
      logical:: test_init
      integer(HID_T):: fid
      integer(HSIZE_T):: global_dims(2)
      integer(HSSIZE_T):: displ(2)

      ! how many particles to generate per process
      nvirt_loc = ceiling(dble(nvirt)/numImages)
      if (thisImage .eq. numImages) nvirt_loc = nvirt - (numImages - 1)*nvirt_loc

      ! Exchanging how many particles each process will read data for
      allocate(nvirt_glob(numImages)[*])
      do n = 0, numImages-1
         otherImage = mod(thisImage+n, numImages)
         if (otherImage==0) otherImage = numImages
         nvirt_glob(thisImage) [otherImage] = nvirt_loc
      end do

      call MPI_INITIALIZED(test_init, ierr)
      if (.not. test_init) call MPI_INIT(ierr)

      ! initializing hdf5
      call h5open_f(ierr)

      call hdf5_parallel_fileopen_read(input_file, fid)

      sync all

      global_dims(1) = dim; global_dims(2) = ntotal
      displ(1) = 0; displ(2) = SUM(nvirt_glob(1:thisImage-1))

      if (nvirt /= 0) call read_particle_data_parallel(thisImage, fid, 'virt', global_dims, displ, &
         parts(ntotal_loc+1:ntotal_loc+nvirt_loc))

      do i = ntotal_loc+1, ntotal_loc+nvirt_loc
         parts(i)%indloc = i
      end do

      ! Closing output file
      call h5fclose_f(fid, ierr)

      ! closing hdf5
      call h5close_f(ierr)

      if (.not. test_init) call MPI_FINALIZE(ierr)

   end subroutine generate_virt_part

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
   subroutine read_particle_data_parallel(thisImage, fid_in, particletype, gdims, ldispl, parts)

      implicit none
      integer(HID_T), intent(in):: fid_in
      character(*), intent(in):: particletype
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
      call hdf5_parallel_read(fid_in, particletype//'/ind', ldispl(2), gdims(2), int_tmp(:, 1))
      parts(:)%indglob = int_tmp(:, 1)
      
      call hdf5_parallel_read(fid_in, particletype//'/type', ldispl(2), gdims(2), int_tmp(:, 1))
      parts(:)%itype = int_tmp(:, 1)

      call hdf5_parallel_read(fid_in, particletype//'/rho', ldispl(2), gdims(2), dbl_tmp(:, 1))
      parts(:)%rho = dbl_tmp(:, 1)

      call hdf5_parallel_read(fid_in, particletype//'/p', ldispl(2), gdims(2), dbl_tmp(:, 1))
      parts(:)%p = dbl_tmp(:, 1)

      deallocate(int_tmp, dbl_tmp)

      allocate( dbl_tmp(dim,nelem))

      write(*,*) thisImage, ldispl
      
      call hdf5_parallel_read(fid_in, particletype//'/x', ldispl, gdims, dbl_tmp)
      do i = 1, nelem
         parts(i)%x(:) = dbl_tmp(:, i)
      end do

      call hdf5_parallel_read(fid_in, particletype//'/v', ldispl, gdims, dbl_tmp)
      do i = 1, nelem
         parts(i)%vx(:) = dbl_tmp(:, i)
      end do


   end subroutine read_particle_data_parallel

end module input_m
