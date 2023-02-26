module ORB_sr_m

   use datatypes, only: particles, system_clock_timer
   use iso_fortran_env, only: lock_type, event_type
   use param_para, only: neighbour_data
   use param, only: tf, f, dim, hsml, halotype

   private
   public:: ORB_sendrecv_diffuse, ORB_sendrecv_halo

contains

   !===================================================================================================================
   subroutine ORB_sendrecv_diffuse(itimestep, thisImage, bounds_loc, repartition_mode, n_process_neighbour, &
                                   neighbours, old_ntotal_loc, ntotal_loc, parts)
      ! Recursive function to exchange physical particles. In cases were subdomain boundaries are updated, the possibility of needing
      ! diffusion is considered

      implicit none
      integer, intent(in):: itimestep, thisImage, repartition_mode, n_process_neighbour
      real(f), intent(in):: bounds_loc(2*dim)
      integer, codimension[*], intent(inout):: ntotal_loc
      type(particles), codimension[*], intent(inout):: parts(:)
      type(neighbour_data), intent(inout):: neighbours(:)
      integer, intent(out):: old_ntotal_loc
      integer:: i, j, k, d, n, searchrange(2), entrydepth, nphys_send_all, ndiffuse, diff_dest, displ0, displ1, &
                neighbourImageIDs(n_process_neighbour)
      real(f):: xmin_loc(dim), xmax_loc(dim), xi(dim), xmin_rem(dim), xmax_rem(dim), dr_min, dr
      real(tf):: tmptime
      logical:: diffuse
      integer, allocatable:: removal_list(:)
      type(lock_type), codimension[*], save:: lock
      type(event_type), codimension[*], save:: neighbourEvent(2)

      ! Initialization
      diffuse = .true.
      entrydepth = 0
      searchrange(:) = [1, ntotal_loc]
      xmin_loc(:) = bounds_loc(1:dim)
      xmax_loc(:) = bounds_loc(dim + 1:2*dim)

      diffuseloop: do while (diffuse)

         ! Searching particles to remove within indices of searchrange(1) and searchrange(2), inclusive.
         ! At node 0, searchrange(1:2) = [1,ntotal_loc)
         allocate (removal_list(searchrange(2) - searchrange(1) + 2))
         nphys_send_all = 0
         do i = searchrange(1), searchrange(2)
            xi(:) = parts(i)%x(:)
            if (any(xi(:) < xmin_loc(:)) .or. any(xi(:) >= xmax_loc(:))) then
               nphys_send_all = nphys_send_all + 1
               removal_list(nphys_send_all) = i
            end if
         end do
         removal_list(nphys_send_all + 1) = 0 ! This is needed due to a quirk in the loop below

         ! If there are any particles that do not belong to the host process,
         ! begin searching for neighbouring processes to send the particle to.
         ! If particle is not contained with subdomain boundaries, send the particle to nearest process neighbouring current host.
         do i = 1, n_process_neighbour
            neighbours(i)%nphys_send = 0
            allocate (neighbours(i)%PhysPackSend(nphys_send_all))
         end do

         tmptime = -system_clock_timer()

         ndiffuse = 0
         if (nphys_send_all > 0) then
            loop_through_parts: do j = 1, nphys_send_all
               i = removal_list(j)
               xi(:) = parts(i)%x(:)
               do n = 1, n_process_neighbour
                  xmin_rem(:) = neighbours(n)%bounds(1:dim)
                  xmax_rem(:) = neighbours(n)%bounds(dim + 1:2*dim)
                  if (.not. (any(xi(:) < xmin_rem(:)) .or. any(xi(:) >= xmax_rem(:)))) then
                     neighbours(n)%nphys_send = neighbours(n)%nphys_send + 1
                     neighbours(n)%PhysPackSend(neighbours(n)%nphys_send) = parts(i)
                     cycle loop_through_parts
                  end if
               end do

               ! if particle makes it past the do loop above, then it belongs to an image that is not the current
               ! image's neighbour
               ndiffuse = ndiffuse + 1
               dr_min = huge(1._f)
               diff_dest = 0
               do n = 1, n_process_neighbour
                  dr = 0._f
                  do d = 1, dim
                     dr = dr + MAX(0._f, &
                                   neighbours(n)%bounds(d) - parts(i)%x(d), &
                                   parts(i)%x(d) - neighbours(n)%bounds(dim + d))**2
                  end do
                  dr = sqrt(dr)
                  if (dr < dr_min) then
                     diff_dest = n
                     dr_min = dr
                  end if
               end do

               neighbours(diff_dest)%nphys_send = neighbours(diff_dest)%nphys_send + 1
               neighbours(diff_dest)%PhysPackSend(neighbours(diff_dest)%nphys_send) = parts(i)

            end do loop_through_parts

         end if

         ! Shifting information to remove sent particles
         n = 0
         if (nphys_send_all > 0 .and. nphys_send_all < ntotal_loc) then
            do i = removal_list(1), ntotal_loc
               if (i == removal_list(n + 1)) then
                  n = n + 1
               else
                  parts(i - n) = parts(i)
               end if
            end do
         end if

         ntotal_loc = ntotal_loc - nphys_send_all
         old_ntotal_loc = ntotal_loc

         ! neighbourImageIDs introduced because sync images didn't like non-contiguous arrays
         ! neighbourImageIDs = neighbours(1:n_process_neighbour)%image
         ! ensuring neighbour images have up-to-date ntotal_loc values
         ! sync all !images(neighbourImageIDs)
         do n = 1, n_process_neighbour
            event post (neighbourEvent(1)[neighbours(n)%image])
         end do
         event wait (neighbourEvent(1), until_count=n_process_neighbour)

         ! each process places particles destined for another image, into that image's particle array
         do i = 1, n_process_neighbour
            if (neighbours(i)%nphys_send > 0) then

               ! lock neighbour image, so this image can update ntotal_loc value
               lock(lock[neighbours(i)%image])
               neighbours(i)%physdispl(1) = ntotal_loc[neighbours(i)%image] + 1
               neighbours(i)%physdispl(2) = neighbours(i)%physdispl(1) + neighbours(i)%nphys_send - 1
               ntotal_loc[neighbours(i)%image] = neighbours(i)%physdispl(2)
               unlock(lock[neighbours(i)%image])
            end if
         end do

         do i = 1, n_process_neighbour
            if (neighbours(i)%nphys_send > 0) then
               parts(neighbours(i)%physdispl(1):neighbours(i)%physdispl(2)) [neighbours(i)%image] = &
                  neighbours(i)%PhysPackSend(1:neighbours(i)%nphys_send)
            end if
         end do

         ! if subdomain boundary update has occurred, check if diffusion is necessary.
         ! Perform if necessary
         if (repartition_mode >= 2) then

            if (thisImage == 1) write (*, '(A)') 'Checking whether diffusion is needed...'

            call co_sum(ndiffuse)

            if (ndiffuse == 0) then
               diffuse = .false.
               if (thisImage == 1) write (*, '(A)') 'No diffusion needed. Continuing...'
            else
               if (thisImage .eq. 1) then
                  write (*, '(A,I0)') 'Diffusion occuring... Current timestep:         ', itimestep
                  write (*, '(A,I0)') '                      Current depth:            ', entrydepth
                  write (*, '(A,I0)') '                      Particles to be diffused: ', ndiffuse
               end if

               entrydepth = entrydepth + 1

               searchrange(:) = [old_ntotal_loc + 1, ntotal_loc]

            end if
         else
            diffuse = .false.
         end if

         deallocate (removal_list)

         ! sync all !images(neighbourImageIDs) ! ensuring transfer between neighbours have completed
         do n = 1, n_process_neighbour
            event post (neighbourEvent(2)[neighbours(n)%image])
         end do
         event wait (neighbourEvent(2), until_count=n_process_neighbour)

         do i = 1, n_process_neighbour
            deallocate (neighbours(i)%PhysPackSend)
         end do

      end do diffuseloop

   end subroutine ORB_sendrecv_diffuse

   !====================================================================================================================
   subroutine ORB_sendrecv_halo(thisImage, bounds_loc, scale_k, n_process_neighbour, neighbours, old_ntotal_loc, &
                                ntotal_loc, nhalo_loc, parts)

      !subroutine responsible for sending sending halo particle information between processes, given
      !predetermined subdomain boundaires.
      !Note: subdomain boundaries are used as inputs (bounds_glob).

      implicit none
      integer, intent(in):: thisImage, n_process_neighbour, old_ntotal_loc
      integer, codimension[*], intent(in):: ntotal_loc
      real(f), intent(in):: bounds_loc(2*dim), scale_k
      type(neighbour_data), intent(inout):: neighbours(n_process_neighbour)
      integer, codimension[*], intent(out):: nhalo_loc
      type(particles), codimension[*], intent(inout):: parts(:)
      integer:: i, j, k, d, n, searchrange(2), neighbourImageIDs(n_process_neighbour), displ0, displ1
      real(f):: xmin_loc(dim), xmax_loc(dim), xi(dim), tmptime
      type(lock_type), codimension[*], save:: lock
      type(event_type), codimension[*], save:: neighbourEvent(2)

      ! initialization
      nhalo_loc = 0

      do i = 1, n_process_neighbour
         neighbours(i)%nhalo_send = 0
         if (allocated(neighbours(i)%halo_pindex)) then
            if (ntotal_loc > size(neighbours(i)%halo_pindex)) deallocate (neighbours(i)%halo_pindex) 
         end if
         if (.not. allocated(neighbours(i)%halo_pindex)) allocate (neighbours(i)%halo_pindex(ntotal_loc))
         if (allocated(neighbours(i)%HaloPackSend)) then
            if (ntotal_loc > size(neighbours(i)%HaloPackSend)) deallocate (neighbours(i)%HaloPackSend)
         end if
         if (.not. allocated(neighbours(i)%HaloPackSend)) allocate (neighbours(i)%HaloPackSend(ntotal_loc))
      end do

      neighbourImageIDs(:) = neighbours(1:n_process_neighbour)%image

      ! halo particle send location determination
      ! first loop loops over currently held particles
      searchrange(:) = [1, old_ntotal_loc]
      tmptime = -system_clock_timer()
      ! begin search
      xmin_loc = bounds_loc(1:dim) + 2._f*scale_k*hsml
      xmax_loc = bounds_loc(dim + 1:2*dim) - 2._f*scale_k*hsml
      do k = 1, 1
         !do i = searchrange(1), searchrange(2)
         do i = 1, ntotal_loc
            xi(:) = parts(i)%x(:)
            if (any([xi(:) <= xmin_loc(:), xi(:) >= xmax_loc(:)])) then
               do j = 1, n_process_neighbour
                  if (all([xi(:) >= neighbours(j)%bounds(1:dim) - 2._f*scale_k*hsml, &
                           xi(:) <= neighbours(j)%bounds(dim + 1:2*dim) + 2._f*scale_k*hsml])) then
                     neighbours(j)%nhalo_send = neighbours(j)%nhalo_send + 1
                     neighbours(j)%halo_pindex(neighbours(j)%nhalo_send) = i
                     neighbours(j)%HaloPackSend(neighbours(j)%nhalo_send) = parts(i)
                     neighbours(j)%HaloPackSend(neighbours(j)%nhalo_send)%itype = &
                        neighbours(j)%HaloPackSend(neighbours(j)%nhalo_send)%itype + sign(halotype, parts(i)%itype)
                  end if
               end do
            end if
         end do

         ! waiting for physicla particles to complete exchange if needed
         ! sync all !images(neighbourImageIDs)
         do n = 1, n_process_neighbour
            event post (neighbourEvent(1)[neighbours(n)%image])
         end do
         event wait (neighbourEvent(1), until_count=n_process_neighbour)

         searchrange(:) = [old_ntotal_loc + 1, ntotal_loc]

      end do

      ! each process places particles destined for another image, into that image's particle array
      do i = 1, n_process_neighbour
         if (neighbours(i)%nhalo_send > 0) then

            ! lock neighbour image, so this image can update ntotal_loc value
            lock(lock[neighbours(i)%image])
            neighbours(i)%halodispl(1) = nhalo_loc[neighbours(i)%image] + 1
            neighbours(i)%halodispl(2) = neighbours(i)%halodispl(1) + neighbours(i)%nhalo_send - 1
            nhalo_loc[neighbours(i)%image] = neighbours(i)%halodispl(2)
            unlock(lock[neighbours(i)%image])
            neighbours(i)%halodispl(:) = neighbours(i)%halodispl(:) + ntotal_loc[neighbours(i)%image]
         end if
      end do

      do i = 1, n_process_neighbour
         if (neighbours(i)%nhalo_send > 0) then
            parts(neighbours(i)%halodispl(1):neighbours(i)%halodispl(2)) [neighbours(i)%image] = &
               neighbours(i)%HaloPackSend(1:neighbours(i)%nhalo_send)
         end if
      end do

      do n = 1, n_process_neighbour
         event post (neighbourEvent(2)[neighbours(n)%image])
      end do
      event wait (neighbourEvent(2), until_count=n_process_neighbour)

      tmptime = tmptime + system_clock_timer()

   end subroutine ORB_sendrecv_halo

end module ORB_sr_m
