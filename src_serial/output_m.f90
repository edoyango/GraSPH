module output_m

   use datatypes, only: particles
   use param, only: output_directory, f, dim

   use hdf5
   use h5lt

   integer(HID_T), private:: file_id, real_group_id, virt_group_id, ghos_group_id
   integer, private:: ierr

contains

   !==============================================================================================================================
   subroutine output(itimestep, save_step, ntotal, nvirt, nghos, parts)
      ! Subroutine to write data to disk. Called intermittently, determined by save_step supplied at run-time

      implicit none
      integer, intent(in):: itimestep, save_step, ntotal, nvirt, nghos
      type(particles), intent(in):: parts(:)
      character(len=4)::number

      write (number, '(I4.4)') itimestep/save_step

      ! Initializing hdf5 interface
      call h5open_f(ierr)

      ! Creating output file
      call h5fcreate_f(trim(output_directory)//"/sph_out"//number//".h5", h5F_ACC_TRUNC_F, file_id, ierr)

      ! Creating groups for each of real, virtual, and ghost particles
      call h5gcreate_f(file_id, "real", real_group_id, ierr)
      call h5gclose_f(real_group_id, ierr)

      call h5gcreate_f(file_id, "virt", virt_group_id, ierr)
      call h5gclose_f(virt_group_id, ierr)

      call h5gcreate_f(file_id, "ghos", ghos_group_id, ierr)
      call h5gclose_f(ghos_group_id, ierr)

      call write_particle_data(file_id, 'real', parts(1:ntotal))
      call write_particle_data(file_id, 'virt', parts(ntotal + 1:ntotal + nvirt))
      call write_particle_data(file_id, 'ghos', parts(ntotal + nvirt + 1:ntotal + nvirt + nghos))

      call h5fclose_f(file_id, ierr)

   end subroutine output

   !===============================================================================================================================
   subroutine write_particle_data(fid_in, group_label, pts_in)

      implicit none
      integer(HID_T), intent(in):: fid_in
      character(*), intent(in):: group_label
      type(particles), intent(in):: pts_in(:)
      integer(HSIZE_T):: data_dims(2)
      integer:: nelem, i
      integer, allocatable:: idatatmp(:, :)
      real(f), allocatable:: fdatatmp(:, :)

      nelem = size(pts_in)

      ! Int data
      data_dims(:) = [nelem, 1]
      allocate (idatatmp(data_dims(1), data_dims(2)))

      ! particle index
      idatatmp(:, 1) = pts_in(:)%ind
      call H5LTmake_dataset_int_f(fid_in, group_label//'/ind', 1, data_dims, idatatmp, ierr)

      ! process ID
      idatatmp(:, 1) = 0
      call H5LTmake_dataset_int_f(fid_in, group_label//'/procid', 1, data_dims, idatatmp, ierr)

      ! particle type
      idatatmp(:, 1) = pts_in(:)%itype
      call H5LTmake_dataset_int_f(fid_in, group_label//'/type', 1, data_dims, idatatmp, ierr)
      deallocate (idatatmp)

      ! float data
      data_dims(:) = [dim, nelem]
      allocate (fdatatmp(data_dims(1), data_dims(2)))

      ! position
      do i = 1, nelem
         fdatatmp(:, i) = pts_in(i)%x(:)
      end do
      select case (f)
      case (kind(1.))
         call H5LTmake_dataset_float_f(fid_in, group_label//'/x', 2, data_dims, fdatatmp, ierr)
      case default
         call H5LTmake_dataset_double_f(fid_in, group_label//'/x', 2, data_dims, fdatatmp, ierr)
      end select

      ! velocity
      do i = 1, nelem
         fdatatmp(:, i) = pts_in(i)%vx(:)
      end do
      select case (f)
      case (kind(1.))
         call H5LTmake_dataset_float_f(fid_in, group_label//'/v', 2, data_dims, fdatatmp, ierr)
      case default
         call H5LTmake_dataset_double_f(fid_in, group_label//'/v', 2, data_dims, fdatatmp, ierr)
      end select

      deallocate (fdatatmp)

      ! density
      data_dims(:) = [nelem, 1]
      allocate (fdatatmp(data_dims(1), data_dims(2)))
      fdatatmp(:, 1) = pts_in(:)%rho
      select case (f)
      case (kind(1.))
         call H5LTmake_dataset_float_f(fid_in, group_label//'/rho', 1, data_dims, fdatatmp, ierr)
      case default
         call H5LTmake_dataset_double_f(fid_in, group_label//'/rho', 1, data_dims, fdatatmp, ierr)
      end select

      ! pressure
      fdatatmp(:, 1) = pts_in(:)%p
      select case (f)
      case (kind(1.))
         call H5LTmake_dataset_float_f(fid_in, group_label//'/p', 1, data_dims, fdatatmp, ierr)
      case default
         call H5LTmake_dataset_double_f(fid_in, group_label//'/p', 1, data_dims, fdatatmp, ierr)
      end select

      deallocate (fdatatmp)

   end subroutine write_particle_data

end module output_m
