!> @file monaghan1994.f90
!> @brief Module and program containing setup code to replicate the classic dambreak experiment from
!>        Monaghan (1994) (https://doi.org/10.1006/jcph.1994.1034). The only known difference is that the Leap-Frog time-integration
!>        is used here instead of the predictor-corrector scheme used in the paper.
!> @author Edward Yang
!> @date 2025/9/22

program main

    use grasph_constants, only: fp, pi
    use grasph_particles, only: particles_container, bp => base_particles
    use weakly_compressible_particles, only: wcp => tait_eos_particles
    use weakly_compressible_interactions, only: fluid_sweeper
    use grasph_pair_sets, only: particle_interactions
    use grasph_time_integration, only: leap_frog_time_integration
    use grasph_kernels, only: grasph_cubic_bspline_kernel
    use grasph_particle_shifting, only: xsph_shifter

    implicit none
    ! parameters to describe geometry
    real(fp), parameter:: dx = 0.04_fp, g = 0._fp, rho0 = 1000._fp
    ! no. of particles in x, y direction in initial geometry of fluid
    integer, parameter:: nfx = 2._fp/dx, nfy = 2._fp/dx
    type(particles_container):: ps(1)
    type(particle_interactions):: pic(1)
    type(grasph_cubic_bspline_kernel):: kernel
    integer:: i, j, k
    type(fluid_sweeper):: sweeper
    type(xsph_shifter):: shifter
    real(fp):: x, y

    ! declare fluid particles
    allocate (wcp::ps(1)%p)

    ! init fluid particles
    select type (ps => ps(1)%p) ! specialise for weakly compressible particles
    class is (wcp)
        call ps%init(n=1976, name="fluid", rho_ref=rho0)
        call ps%register_x%register(ps%ps(1), ps%ps(1)%x, ps%ps(1)%v)
        call ps%register_v%register(ps%ps(1), ps%ps(1)%v, ps%ps(1)%dvxdt)
        call ps%register_v%register(ps%ps(1), ps%ps(1)%rho, ps%ps(1)%drhodt)
    end select

    ! initialize geometry
    ! generate particles in a grid and save only the ones within the 1 radius circle
    k = 0
    do i = 0, nfx - 1
        do j = 0, nfy - 1
            x = -1._fp + (i + 0.5_fp)*dx
            y = -1._fp + (j + 0.5_fp)*dx
            if (x*x + y*y < 1._fp) then
                k = k + 1
                ps(1)%p%ps(k)%id = k
                ps(1)%p%ps(k)%type = 1
                ps(1)%p%ps(k)%x(1) = x
                ps(1)%p%ps(k)%x(2) = y
                ps(1)%p%ps(k)%c = 1400._fp
                ps(1)%p%ps(k)%rho = rho0
                ps(1)%p%ps(k)%mass = pi*rho0/1976
                ps(1)%p%ps(k)%v(1) = -100._fp*x
                ps(1)%p%ps(k)%v(2) = 100._fp*y
            end if
        end do
    end do

    ! init interactions
    sweeper%artvisc_alpha = 0.01_fp
    sweeper%artvisc_beta = 0._fp
    sweeper%h = 1.2_fp*dx
    sweeper%g = g
    shifter%epsilon = 0.5_fp
    shifter%update_rhs = .true.
    call pic(1)%init(30, ps(1)%p, sweeper=sweeper, shifter=shifter)

    ! init kernel
    call kernel%init(2, 1.2_fp*dx)

    call leap_frog_time_integration(5000, 10, 10, ps, pic, 0.05_fp, kernel, "/home/edwardy/test", output_comp_level=4)

end program main
