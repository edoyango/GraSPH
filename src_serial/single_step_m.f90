module single_step_m

   public:: single_step

contains

   !==============================================================================================================================
   pure subroutine single_step(ki, ntotal, nvirt, nghos, niac, pairs, parts, dvxdti, drhoi, nexti)
      ! Container subroutine for all the rate-of-change calculations. Rate-of-changes are calculated seperately and then summed as
      ! required

      use datatypes, only: particles, interactions
      use param, only: dim, f, g

      use input_m, only: virt_mirror
      use material_rates_m, only: art_visc, int_force, con_density, art_visc_coeff, int_force_coeff

      implicit none
      integer, intent(in):: ki, ntotal, nvirt, nghos, niac, nexti(:)
      type(interactions), intent(in):: pairs(:)
      type(particles), intent(inout):: parts(:)
      real(f), intent(out):: dvxdti(dim, ntotal), drhoi(ntotal)
      integer:: i, j, k
      real(f):: a_coeff
      type(particles):: p_i

      drhoi(1:ntotal) = 0._f
      do i = 1,ntotal
         dvxdti(1:dim-1,i) = 0._f
         dvxdti(dim,i) = -g
      end do

      ! looping through interaction pairs to calculate forces/density change
      do i = 1,ntotal+nvirt+nghos
         
         do k = nexti(i),nexti(i+1)-1
         
            j = pairs(k)%j
            
            if (parts(i)%itype > 0 .and. parts(j)%itype < 0) then
               call virt_mirror(parts(i), parts(j))
            elseif (parts(i)%itype < 0 .and. parts(j)%itype > 0) then
               call virt_mirror(parts(j), parts(i))
            end if
            
            !Density approximation or change rate
            call con_density(ki, parts(i), parts(j), pairs(k)%dwdx, drhoi(i), drhoi(j))
            
            ! calculating coefficients for pressure force and artificial viscosity
            a_coeff = int_force_coeff(ki,parts(i),parts(j)) + art_visc_coeff(ki,parts(i),parts(j)) 
            dvxdti(:,i) = dvxdti(:,i) + pairs(k)%dwdx(:)*a_coeff                                           
            dvxdti(:,j) = dvxdti(:,j) - pairs(k)%dwdx(:)*a_coeff
            
         end do

      end do
            
   end subroutine single_step

end module single_step_m
