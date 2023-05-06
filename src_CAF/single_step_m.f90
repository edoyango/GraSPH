module single_step_m

   use datatypes, only: particles, interactions, system_clock_timer
   use kernel_m, only: kernel_dwdx
   use param, only: dim, f, g, hsml
   use material_rates_m, only: con_density, art_visc_coeff, int_force_coeff

   private
   public:: single_step

contains

   !==============================================================================================================================
   pure subroutine single_step(parts, niac, pairs, dvxdt, drhodt, nexti)
      ! Container subroutine for all the rate-of-change calculations. Rate-of-changes are calculated seperately and then summed as
      ! required

      !    use ORB_sr_m, only: ORB_sendrecv_haloupdate

      implicit none
      integer, intent(in):: niac, nexti(:)
      type(particles), intent(inout):: parts
      type(interactions), intent(in):: pairs(:)
      real(f), intent(out):: dvxdt(:, :), drhodt(:)
      integer:: i, j, k
      real(f):: a_coeff, tdwdx(dim)
      real(f), allocatable:: prho(:)

      allocate (prho(parts%nsum()))

      drhodt(1:parts%ntotal_loc+parts%nvirt_loc) = 0._f
      do i = 1, parts%ntotal_loc+parts%nvirt_loc
         dvxdt(1:dim - 1, i) = 0._f
         dvxdt(dim, i) = -g
      end do

      do i = 1, parts%nsum()
         prho(i) = parts%p(i)/parts%rho(i)**2
      end do

      do i = 1, parts%nsum()
         do k = nexti(i), nexti(i + 1) - 1

            j = pairs(k)%j

            tdwdx(:) = kernel_dwdx(pairs(k)%dx, hsml)

            !Density approximation or change rate
            call con_density(parts%vx(:, i), parts%vx(:, j), tdwdx, drhodt(i), drhodt(j))

            ! calculating coefficients for pressure force and artificial viscosity
            a_coeff = int_force_coeff(prho(i), prho(j)) + &
                      art_visc_coeff(parts%rho(i), parts%vx(:, i), parts%rho(j), parts%vx(:, j), pairs(k)%dx)
            dvxdt(:, i) = dvxdt(:, i) + tdwdx(:)*a_coeff
            dvxdt(:, j) = dvxdt(:, j) - tdwdx(:)*a_coeff

         end do

      end do

   end subroutine single_step

end module single_step_m
