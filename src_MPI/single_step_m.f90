module single_step_m

   use datatypes, only: particles, interactions

   private
   public:: single_step

contains

   !==============================================================================================================================
   pure subroutine single_step(ntotal_loc, nhalo_loc, nvirt_loc, parts, niac, pairs, dvxdti, drhoi, nexti)
      ! Container subroutine for all the rate-of-change calculations. Rate-of-changes are calculated seperately and then summed as
      ! required

      use mpi_f08, only: MPI_WTIME
      use param, only: dim, f, g

      use material_rates_m, only: int_force, art_visc, con_density, ext_force, art_visc_coeff, int_force_coeff
      use ORB_sr_m, only: ORB_sendrecv_haloupdate

      implicit none
      integer, intent(in):: niac, ntotal_loc, nhalo_loc, nvirt_loc, nexti(:)
      type(particles), intent(inout):: parts(:)
      type(interactions), intent(in):: pairs(:)
      real(f), intent(out):: dvxdti(:, :), drhoi(:)
      integer:: i, j, k
      real(f):: a_coeff
      real(f), allocatable:: prho(:)

      allocate (prho(ntotal_loc + nhalo_loc + nvirt_loc))

      drhoi(1:ntotal_loc) = 0._f
      do i = 1, ntotal_loc
         dvxdti(1:dim - 1, i) = 0._f
         dvxdti(dim, i) = -g
      end do

      do i = 1, ntotal_loc + nhalo_loc + nvirt_loc
         prho(i) = parts(i)%p/parts(i)%rho**2
      end do

      do i = 1, ntotal_loc + nhalo_loc + nvirt_loc
         do k = nexti(i), nexti(i + 1) - 1

            j = pairs(k)%j

            !Density approximation or change rate
            call con_density(parts(i), parts(j), pairs(k)%dwdx, drhoi(i), drhoi(j))

            ! calculating coefficients for pressure force and artificial viscosity
            a_coeff = int_force_coeff(prho(i), prho(j)) + art_visc_coeff(parts(i), parts(j))
            dvxdti(:, i) = dvxdti(:, i) + pairs(k)%dwdx(:)*a_coeff
            dvxdti(:, j) = dvxdti(:, j) - pairs(k)%dwdx(:)*a_coeff

         end do

      end do

   end subroutine single_step

end module single_step_m
