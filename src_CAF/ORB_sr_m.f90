module ORB_sr_m

   use datatypes, only: particles
   use iso_fortran_env, only: lock_type
   use param_para, only: neighbour_data
   use param, only: f, dim, hsml

   private
   public:: ORB_sendrecv_diffuse

contains

   !===================================================================================================================
   subroutine ORB_sendrecv_diffuse(itimestep, thisImage, bounds_loc, repartition_mode, n_process_neighbour, &
                                    neighbours, n_recv_all, ntotal_loc, parts)
      ! Recursive function to exchange physical particles. In cases were subdomain boundaries are updated, the possibility of needing
      ! diffusion is considered

      implicit none
      integer, intent(in):: itimestep, thisImage, repartition_mode, n_process_neighbour
      real(f), intent(in):: bounds_loc(2*dim)
      integer, codimension[*], intent(inout):: ntotal_loc
      type(particles), codimension[*], intent(inout):: parts(:)
      type(neighbour_data), intent(inout):: neighbours(:)
      integer, intent(out):: n_recv_all
      integer:: i, j, k, d, n, searchrange(2), entrydepth=0, nphys_send_all, ndiffuse, diff_dest, displ0, displ1, &
         old_ntotal_loc, neighbourImageIDs(n_process_neighbour)
      real(f):: xmin_loc(dim), xmax_loc(dim), xi(dim), xmin_rem(dim), xmax_rem(dim), dr_min, dr
      logical:: diffuse = .true.
      integer, allocatable:: removal_list(:)
      type(lock_type), save:: lock[*]

      ! Initialization
      searchrange(:) = [1, ntotal_loc]
      xmin_loc(:) = bounds_loc(1:dim)
      xmax_loc(:) = bounds_loc(dim+1:2*dim)

      diffuseloop: do while (diffuse)

         ! Searching particles to remove within indices of searchrange(1) and searchrange(2), inclusive.
         ! At node 0, searchrange(1:2) = [1,ntotal_loc)
         allocate( removal_list(searchrange(2)-searchrange(1)+2))
         nphys_send_all = 0
         do i = searchrange(1), searchrange(2)
            xi(:) = parts(i)%x(:)
            if (any([xi(:) < xmin_loc(:), xi(:) >= xmax_loc(:)])) then
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
            allocate( neighbours(i)%PhysPackSend(nphys_send_all))
         end do

         ndiffuse = 0
         if (nphys_send_all > 0) then
            loop_through_parts: do j = 1, nphys_send_all
               i = removal_list(j)
               xi(:) = parts(i)%x(:)
               do n = 1, n_process_neighbour
                  xmin_rem(:) = neighbours(n)%bounds(1:dim)
                  xmax_rem(:) = neighbours(n)%bounds(dim+1:2*dim)
                  if (.not. any([xi(:) < xmin_rem(:), xi(:) >= xmax_rem(:)])) then
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
                  do d = 1,dim
                     dr = dr + MAX(0._f, &
                                   neighbours(n)%bounds(d) - parts(i)%x(d), &
                                   parts(i)%x(d) - neighbours(n)%bounds(dim+d))**2
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
         if ( nphys_send_all > 0 .and. nphys_send_all < ntotal_loc) then
            do i = removal_list(1), ntotal_loc
               if (i == removal_list(n+1)) then
                  n = n + 1
               else
                  parts(i - n) = parts(i)
               end if
            end do
         end if

         ntotal_loc = ntotal_loc - nphys_send_all
         old_ntotal_loc = ntotal_loc

         ! neighbourImageIDs introduced because sync images didn't like non-contiguous arrays
         neighbourImageIDs = neighbours(1:n_process_neighbour)%image
         ! ensuring neighbour images have up-to-date ntotal_loc values
         sync images(neighbourImageIDs)

         ! each process places particles destined for another image, into that image's particle array
         do i = 1, n_process_neighbour
            if (neighbours(i)%nphys_send > 0) then

               ! lock neighbour image, so this image can update ntotal_loc value
               lock(lock[neighbours(i)%image])
               displ0 = ntotal_loc[neighbours(i)%image] + 1
               displ1 = displ0 + neighbours(i)%nphys_send - 1
               ntotal_loc[neighbours(i)%image] = displ1
               unlock(lock[neighbours(i)%image])

               parts(displ0:displ1)[neighbours(i)%image] = neighbours(i)%PhysPackSend(1:neighbours(i)%nphys_send)
            end if
         end do

         ! if subdomain boundary update has occurred, check if diffusion is necessary.
         ! Perform if necessary
         if (repartition_mode >= 2) then

            if (thisImage == 1) write(*, '(A)') 'Checking whether diffusion is needed...'

            call co_sum(ndiffuse)

            if (ndiffuse == 0) then
               diffuse = .false.
               if (thisImage==1) write(*, '(A)') 'No diffusion needed. Continuing...'
            else
               if (thisImage .eq. 1) then
                  write (*, '(A,I0)') 'Diffusion occuring... Current timestep:         ', itimestep
                  write (*, '(A,I0)') '                      Current depth:            ', entrydepth
                  write (*, '(A,I0)') '                      Particles to be diffused: ', ndiffuse
               end if

               entrydepth = entrydepth + 1

               searchrange(:) = [old_ntotal_loc+1, ntotal_loc]
               
            end if
         else
            diffuse = .false.
         end if

         deallocate(removal_list)

         do i = 1, n_process_neighbour
            deallocate( neighbours(i)%PhysPackSend)
         end do

         sync images(neighbourImageIDs) ! ensuring transfer between neighbours have completed
            
      end do diffuseloop

   end subroutine ORB_sendrecv_diffuse

end module ORB_sr_m