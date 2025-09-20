!> @file grasph_pairs.f90
!> @brief Module containing types and subroutines related to finding interacting pairs of particles.
!> @author Edward Yang
!> @date 2025-06-09
module grasph_pairs

    use grasph_constants, only: fp
    use grasph_kernels, only: grasph_base_kernel

    implicit none

    private

    !> @brief Type to hold information about sets of pairs of particles.
    type particle_pairs
        !> @brief The number of LHS particles.
        integer:: n = 0
        !> @brief Maximum number of particle interactions expected for each LHS particle.
        integer:: npairs_per_particle = 0
        !> @brief Total number of particle pairs found.
        integer:: npairs_total = 0
        !> @brief Number of spatial dimensions.
        integer:: ndims = 0
        !> @brief ij particle pair indices.
        integer, allocatable:: pair_ij(:, :)
        real(fp), allocatable:: w(:)
        !> @brief The kernel gradients calculated for the given pair.
        real(fp), allocatable:: dwdx(:, :)
        !> @brief Whether the particle_pairs has been initialized.
        logical:: initialized = .false.
    contains
        procedure:: init => particle_pairs_init
    end type particle_pairs

    ! ! abstract interface provided to be used in function pointers
    ! abstract interface
    !     pure subroutine find_pairs_fixed_h(x, ndims, n, cutoff, kernel, npairs_per_particle, pairs)
    !         import:: fp, grasph_base_kernel, particle_pairs
    !         type(particle_pairs), intent(inout):: pairs
    !         integer, intent(in):: ndims, n, npairs_per_particle
    !         real(fp), intent(in):: x(ndims, n), cutoff
    !         class(grasph_base_kernel), intent(in):: kernel
    !     end subroutine find_pairs_fixed_h
    ! end interface

    !> @brief The overloaded interface to find pairs of particles using the direct search (aka
    !>        brute-force) algorithm. Overloaded with a version that finds pairs amongst a single
    !>        set of particles, and two sets of particles. In both instances, particles are assumed
    !>        to have the same cutoff distance.
    interface dsearch
        module procedure dsearch_self, dsearch_other
    end interface dsearch

    !> @brief The overloaded interface to find pairs of particles using the cell-lists algorithm.
    !>        Overloaded with a version that finds pairs amongst a single set of particles, and two
    !>        sets of particles. In both instances, particles are assumed to have the same cutoff
    !>        distance.
    interface cell_list_search
        module procedure cell_list_search_self, cell_list_search_other
    end interface cell_list_search

    public:: particle_pairs, dsearch, cell_list_search

contains

    !> @brief Subroutine to initiliaze the particle_pairs derived type.
    !> @param self The particle_pairs instance to initialize.
    !> @param n The number of LHS particles involved in the search.
    !> @param npairs_per_particle The maximum number of pairs each particle will have.
    !> @param ndims The number of dimensions of the problem.
    pure subroutine particle_pairs_init(self, n, npairs_per_particle, ndims)

        class(particle_pairs), intent(inout):: self
        integer, intent(in):: n, npairs_per_particle, ndims

        self%ndims = ndims
        self%n = n
        self%npairs_per_particle = npairs_per_particle
        self%npairs_total = 0

        if (self%initialized) deallocate (self%pair_ij, self%w, self%dwdx)
        allocate (self%pair_ij(2, n*npairs_per_particle))
        allocate (self%w(n*npairs_per_particle), self%dwdx(ndims, n*npairs_per_particle))

        self%initialized = .true.

    end subroutine particle_pairs_init

    !> @brief The direct-search fixed-radius neighbour search algorithm. Finds pairs within a single
    !>        set of particles.
    !> @param x The positions of the particles.
    !> @param cutoff The cutoff distance to find pairs of particles within.
    !> @param kernel The SPH kernel to calculate values and gradient values with.
    !> @param pairs The particle_pairs instance to populate with the search.
    pure subroutine dsearch_self(x, cutoff, kernel, pairs)

        type(particle_pairs), intent(inout):: pairs
        real(fp), intent(in):: x(pairs%ndims, pairs%n), cutoff
        class(grasph_base_kernel), intent(in):: kernel
        integer:: i, j
        real(fp):: dx(pairs%ndims)

        pairs%npairs_total = 0

        do i = 1, pairs%n - 1
            do j = i + 1, pairs%n
                dx(:) = x(:, i) - x(:, j)
                if (sum(dx(:)**2) < cutoff*cutoff) then
                    pairs%npairs_total = pairs%npairs_total + 1
                    pairs%pair_ij(1, pairs%npairs_total) = i
                    pairs%pair_ij(2, pairs%npairs_total) = j
                    call kernel%values(dx, pairs%w(pairs%npairs_total), pairs%dwdx(:, pairs%npairs_total))
                end if
            end do
        end do
    end subroutine dsearch_self

    !> @brief The direct-search fixed-radius neighbour search algorithm. Finds pairs between two
    !>        distinct sets of particles. Kernel values and gradients are always calculated with
    !>        respect to the LHS particles.
    !> @param x_lhs The positions of the LHS particles.
    !> @param x_rhs The positions of the RHS particles.
    !> @param n_rhs The number of RHS particles.
    !> @param cutoff The cutoff distance to find pairs of particles within.
    !> @param kernel The SPH kernel to calculate values and gradient values with.
    !> @param pairs The particle_pairs instance to populate with the search.
    pure subroutine dsearch_other(x_lhs, x_rhs, n_rhs, cutoff, kernel, pairs)

        integer, intent(in):: n_rhs
        type(particle_pairs), intent(inout):: pairs
        real(fp), intent(in):: x_lhs(pairs%ndims, pairs%n), x_rhs(pairs%ndims, n_rhs), cutoff
        class(grasph_base_kernel), intent(in):: kernel
        integer:: i, j
        real(fp):: dx(pairs%ndims)

        pairs%npairs_total = 0

        do i = 1, pairs%n
            do j = 1, n_rhs
                dx(:) = x_lhs(:, i) - x_rhs(:, j)
                if (sum(dx(:)**2) < cutoff*cutoff) then
                    pairs%npairs_total = pairs%npairs_total + 1
                    pairs%pair_ij(1, pairs%npairs_total) = i
                    pairs%pair_ij(2, pairs%npairs_total) = j
                    call kernel%values(dx, pairs%w(pairs%npairs_total), pairs%dwdx(:, pairs%npairs_total))
                end if
            end do
        end do
    end subroutine dsearch_other

    !> @brief The cell-lists fixed-radius neighbour search algorithm. Finds pairs within a single
    !>        set of particles.
    !> @param x The positions of the particles.
    !> @param cutoff The cutoff distance to find pairs of particles within.
    !> @param kernel The SPH kernel to calculate values and gradient values with.
    !> @param pairs The particle_pairs instance to populate with the search.
    pure subroutine cell_list_search_self(x, cutoff, kernel, pairs)

        type(particle_pairs), intent(inout):: pairs
        real(fp), intent(in):: x(pairs%ndims, pairs%n), cutoff
        class(grasph_base_kernel), intent(in):: kernel
        real(fp):: minextents(pairs%ndims), maxextents(pairs%ndims), dcell
        integer:: i, ngridx(pairs%ndims), grid_idx(pairs%ndims, pairs%n) ! might need to be allocatable in the future...

        ! define grid
        dcell = cutoff ! not functionally meaningful, but helpful conceptually
        minextents(:) = minval(x, 2) - 2._fp*dcell
        maxextents(:) = maxval(x, 2) + 2._fp*dcell
        ngridx(:) = int((maxextents(:) - minextents(:))/dcell) + 1
        ! technically, maxextents should be adjusted, but it isn't used from herin

        do i = 1, pairs%n
            grid_idx(:, i) = int((x(:, i) - minextents(:))/dcell) + 1
        end do

        pairs%npairs_total = 0

        select case (pairs%ndims)
        case (2)
            call grid_sweep_2d_self(cutoff, kernel, ngridx, grid_idx, x, pairs)
        case (3)
            call grid_sweep_3d_self(cutoff, kernel, ngridx, grid_idx, x, pairs)
        case default
            error stop "cell_list_search: only 2d and 3d cases are supported!"
        end select

    end subroutine cell_list_search_self

    !> @brief The 2d version of the grid sweep in the cell-lists pair-finding single particle set
    !>        pair search.
    !> @param cutoff The cutoff distance to find pairs of particles within.
    !> @param kernel The SPH kernel to calculate values and gradient values with.
    !> @param ngridx The number of grid cells in each dimension.
    !> @param grid_idx The grid cells that each particle in the set belongs to.
    !> @param x The positions of the particles.
    !> @param pairs The particle_pairs instance to populate with the search.
    pure subroutine grid_sweep_2d_self(cutoff, kernel, ngridx, grid_idx, x, pairs)

        type(particle_pairs), intent(inout):: pairs
        real(fp), intent(in):: cutoff, x(2, pairs%n)
        class(grasph_base_kernel), intent(in):: kernel
        integer, intent(in):: ngridx(2), grid_idx(2, pairs%n)
        integer:: i, j, icell, jcell, jj, pic
        real(fp):: dx(2)
        integer, allocatable:: n_in_cell(:, :), p_in_cell(:, :, :)

        allocate (n_in_cell(ngridx(1), ngridx(2)), source=0)
        allocate (p_in_cell(pairs%npairs_per_particle, ngridx(1), ngridx(2)))

        ! populate grid
        do i = 1, pairs%n
            icell = grid_idx(1, i)
            jcell = grid_idx(2, i)
            n_in_cell(icell, jcell) = n_in_cell(icell, jcell) + 1
            p_in_cell(n_in_cell(icell, jcell), icell, jcell) = i
        end do

        ! sweep half adjacent cells
        do i = 1, pairs%n
            ! current cell
            icell = grid_idx(1, i)
            jcell = grid_idx(2, i)
            ! special sweep because of j>i check
            do pic = 1, n_in_cell(icell, jcell)
                j = p_in_cell(pic, icell, jcell)
                if (j > i) then
                    dx(:) = x(:, i) - x(:, j)
                    if (sum(dx*dx) < cutoff*cutoff) then
                        pairs%npairs_total = pairs%npairs_total + 1
                        pairs%pair_ij(1, pairs%npairs_total) = i
                        pairs%pair_ij(2, pairs%npairs_total) = j
                        call kernel%values(dx, pairs%w(pairs%npairs_total), pairs%dwdx(:, pairs%npairs_total))
                    end if
                end if
            end do
            ! right cell
            icell = icell + 1
            call sweep_cell(cutoff, 2, x(:, i), pairs%n, x, n_in_cell(icell, jcell), &
                            p_in_cell(:, icell, jcell), kernel, pairs, i)
            ! top row
            jcell = jcell + 1
            do icell = grid_idx(1, i) - 1, grid_idx(1, i) + 1
                call sweep_cell(cutoff, 2, x(:, i), pairs%n, x, n_in_cell(icell, jcell), &
                                p_in_cell(:, icell, jcell), kernel, pairs, i)
            end do
        end do

        deallocate (n_in_cell, p_in_cell)

    end subroutine grid_sweep_2d_self

    !> @brief The 3d version of the grid sweep in the cell-lists pair-finding single particle set
    !>        pair search.
    !> @param cutoff The cutoff distance to find pairs of particles within.
    !> @param kernel The SPH kernel to calculate values and gradient values with.
    !> @param ngridx The number of grid cells in each dimension.
    !> @param grid_idx The grid cells that each particle in the set belongs to.
    !> @param x The positions of the particles.
    !> @param pairs The particle_pairs instance to populate with the search.
    pure subroutine grid_sweep_3d_self(cutoff, kernel, ngridx, grid_idx, x, pairs)

        type(particle_pairs), intent(inout):: pairs
        real(fp), intent(in):: cutoff, x(3, pairs%n)
        class(grasph_base_kernel), intent(in):: kernel
        integer, intent(in):: ngridx(3), grid_idx(3, pairs%n)
        integer:: i, j, icell, jcell, kcell, jj, pic
        real(fp):: dx(3)
        integer, allocatable:: n_in_cell(:, :, :), p_in_cell(:, :, :, :)

        allocate (n_in_cell(ngridx(1), ngridx(2), ngridx(3)), source=0)
        allocate (p_in_cell(pairs%npairs_per_particle, ngridx(1), ngridx(2), ngridx(3)))

        ! populate grid
        n_in_cell(:, :, :) = 0
        do i = 1, pairs%n
            icell = grid_idx(1, i)
            jcell = grid_idx(2, i)
            kcell = grid_idx(3, i)
            n_in_cell(icell, jcell, kcell) = n_in_cell(icell, jcell, kcell) + 1
            p_in_cell(n_in_cell(icell, jcell, kcell), icell, jcell, kcell) = i
        end do

        ! sweep half adjacent cells
        do i = 1, pairs%n
            ! current cell
            icell = grid_idx(1, i)
            jcell = grid_idx(2, i)
            kcell = grid_idx(3, i)
            ! special sweep because of j>i check
            do pic = 1, n_in_cell(icell, jcell, kcell)
                j = p_in_cell(pic, icell, jcell, kcell)
                if (j > i) then
                    dx(:) = x(:, i) - x(:, j)
                    if (sum(dx*dx) < cutoff*cutoff) then
                        pairs%npairs_total = pairs%npairs_total + 1
                        pairs%pair_ij(1, pairs%npairs_total) = i
                        pairs%pair_ij(2, pairs%npairs_total) = j
                        call kernel%values(dx, pairs%w(pairs%npairs_total), pairs%dwdx(:, pairs%npairs_total))
                    end if
                end if
            end do
            ! right cell
            icell = icell + 1
            call sweep_cell(cutoff, 3, x(:, i), pairs%n, x, n_in_cell(icell, jcell, kcell), &
                            p_in_cell(:, icell, jcell, kcell), kernel, pairs, i)
            ! north-middle layer
            jcell = jcell + 1
            do icell = grid_idx(1, i) - 1, grid_idx(1, i) + 1
                call sweep_cell(cutoff, 3, x(:, i), pairs%n, x, n_in_cell(icell, jcell, kcell), &
                                p_in_cell(:, icell, jcell, kcell), kernel, pairs, i)
            end do
            ! top layer
            kcell = kcell + 1
            do jcell = grid_idx(2, i) - 1, grid_idx(2, i) + 1
                do icell = grid_idx(1, i) - 1, grid_idx(1, i) + 1
                    call sweep_cell(cutoff, 3, x(:, i), pairs%n, x, n_in_cell(icell, jcell, kcell), &
                                    p_in_cell(:, icell, jcell, kcell), kernel, pairs, i)
                end do
            end do
        end do

        deallocate (n_in_cell, p_in_cell)

    end subroutine grid_sweep_3d_self

    !> @brief The cell-lists fixed-radius neighbour search algorithm. Finds pairs between two
    !>        distinct sets of particles. Kernel values and gradients are always calculated with
    !>        respect to the LHS particles.
    !> @param x_lhs The positions of the LHS particles.
    !> @param x_rhs The positions of the RHS particles.
    !> @param n_rhs The number of RHS particles.
    !> @param cutoff The cutoff distance to find pairs of particles within.
    !> @param kernel The SPH kernel to calculate values and gradient values with.
    !> @param pairs The particle_pairs instance to populate with the search.
    pure subroutine cell_list_search_other(x_lhs, x_rhs, n_rhs, cutoff, kernel, pairs)

        integer, intent(in):: n_rhs
        type(particle_pairs), intent(inout):: pairs
        real(fp), intent(in):: x_lhs(pairs%ndims, pairs%n), x_rhs(pairs%ndims, n_rhs), cutoff
        class(grasph_base_kernel), intent(in):: kernel
        real(fp):: minextents(pairs%ndims), maxextents(pairs%ndims), dcell
        integer:: i, ngridx(pairs%ndims), grid_idx(pairs%ndims, n_rhs) ! might need to be allocatable in the future...

        ! define grid
        dcell = cutoff ! not functionally meaningful, but helpful conceptually
        minextents(:) = min(minval(x_lhs, 2), minval(x_rhs, 2)) - 2._fp*dcell
        maxextents(:) = max(maxval(x_lhs, 2), maxval(x_rhs, 2)) + 2._fp*dcell
        ngridx(:) = int((maxextents(:) - minextents(:))/dcell) + 1
        ! technically, maxextents should be adjusted, but it isn't used from herin

        do i = 1, n_rhs
            grid_idx(:, i) = int((x_rhs(:, i) - minextents(:))/dcell) + 1
        end do

        pairs%npairs_total = 0

        select case (pairs%ndims)
        case (2)
            call grid_sweep_2d_other(cutoff, kernel, minextents, ngridx, grid_idx, x_lhs, x_rhs, n_rhs, pairs)
        case (3)
            call grid_sweep_3d_other(cutoff, kernel, minextents, ngridx, grid_idx, x_lhs, x_rhs, n_rhs, pairs)
        case default
            error stop "cell_list_search: only 2d and 3d cases are supported!"
        end select

    end subroutine cell_list_search_other

    !> @brief The 2d version of the grid sweep in the cell-lists pair-finding two particle set
    !>        pair search.
    !> @param kernel The SPH kernel to calculate values and gradient values with.
    !> @param minextents The minimum coordinates of the grid (corresponding to the
    !>                   left-/south-/bottom-most coordinate)
    !> @param ngridx The number of grid cells in each dimension.
    !> @param grid_idx The grid cells that each particle in the RHS set belongs to.
    !> @param x_lhs The positions of the LHS particles.
    !> @param x_rhs The positions of the RHS particles.
    !> @param n_rhs The number of RHS particles.
    !> @param pairs The particle_pairs instance to populate with the search.
    pure subroutine grid_sweep_2d_other(cutoff, kernel, minextents, ngridx, grid_idx, x_lhs, x_rhs, n_rhs, pairs)

        integer, intent(in):: n_rhs
        type(particle_pairs), intent(inout):: pairs
        real(fp), intent(in):: cutoff, minextents(2), x_lhs(2, pairs%n), x_rhs(2, n_rhs)
        class(grasph_base_kernel), intent(in):: kernel
        integer, intent(in):: ngridx(2), grid_idx(2, n_rhs)
        integer:: i, icell, jcell, this_cell(2)
        real(fp):: dx(2)
        integer, allocatable:: n_in_cell(:, :), p_in_cell(:, :, :)

        allocate (n_in_cell(ngridx(1), ngridx(2)), source=0)
        allocate (p_in_cell(pairs%npairs_per_particle, ngridx(1), ngridx(2)))

        ! populate grid
        do i = 1, n_rhs
            icell = grid_idx(1, i)
            jcell = grid_idx(2, i)
            n_in_cell(icell, jcell) = n_in_cell(icell, jcell) + 1
            p_in_cell(n_in_cell(icell, jcell), icell, jcell) = i
        end do

        ! sweep all adjacent cells
        do i = 1, pairs%n
            this_cell(:) = int((x_lhs(:, i) - minextents(:))/cutoff) + 1
            do jcell = this_cell(2) - 1, this_cell(2) + 1
                do icell = this_cell(1) - 1, this_cell(1) + 1
                    call sweep_cell(cutoff, 2, x_lhs(:, i), n_rhs, x_rhs, n_in_cell(icell, jcell), &
                                    p_in_cell(:, icell, jcell), kernel, pairs, i)
                end do
            end do
        end do

        deallocate (n_in_cell, p_in_cell)

    end subroutine grid_sweep_2d_other

    !> @brief The 3d version of the grid sweep in the cell-lists pair-finding two particle set
    !>        pair search.
    !> @param kernel The SPH kernel to calculate values and gradient values with.
    !> @param minextents The minimum coordinates of the grid (corresponding to the
    !>                   left-/south-/bottom-most coordinate)
    !> @param ngridx The number of grid cells in each dimension.
    !> @param grid_idx The grid cells that each particle in the RHS set belongs to.
    !> @param x_lhs The positions of the LHS particles.
    !> @param x_rhs The positions of the RHS particles.
    !> @param n_rhs The number of RHS particles.
    !> @param pairs The particle_pairs instance to populate with the search.
    pure subroutine grid_sweep_3d_other(cutoff, kernel, minextents, ngridx, grid_idx, x_lhs, x_rhs, n_rhs, pairs)

        integer, intent(in):: n_rhs
        type(particle_pairs), intent(inout):: pairs
        real(fp), intent(in):: minextents(3), cutoff, x_lhs(3, pairs%n), x_rhs(3, n_rhs)
        class(grasph_base_kernel), intent(in):: kernel
        integer, intent(in):: ngridx(3), grid_idx(3, n_rhs)
        integer:: i, icell, jcell, kcell, this_cell(3)
        real(fp):: dx(3)
        integer, allocatable:: n_in_cell(:, :, :), p_in_cell(:, :, :, :)

        allocate (n_in_cell(ngridx(1), ngridx(2), ngridx(3)), source=0)
        allocate (p_in_cell(pairs%npairs_per_particle, ngridx(1), ngridx(2), ngridx(3)))

        ! populate grid
        n_in_cell(:, :, :) = 0
        do i = 1, n_rhs
            icell = grid_idx(1, i)
            jcell = grid_idx(2, i)
            kcell = grid_idx(3, i)
            n_in_cell(icell, jcell, kcell) = n_in_cell(icell, jcell, kcell) + 1
            p_in_cell(n_in_cell(icell, jcell, kcell), icell, jcell, kcell) = i
        end do

        ! sweep all adjacent cells
        do i = 1, pairs%n
            this_cell(:) = int((x_lhs(:, i) - minextents(:))/cutoff) + 1
            do kcell = this_cell(3) - 1, this_cell(3) + 1
                do jcell = this_cell(2) - 1, this_cell(2) + 1
                    do icell = this_cell(1) - 1, this_cell(1) + 1
                        call sweep_cell(cutoff, 3, x_lhs(:, i), n_rhs, x_rhs, n_in_cell(icell, jcell, kcell), &
                                        p_in_cell(:, icell, jcell, kcell), kernel, pairs, i)
                    end do
                end do
            end do
        end do

        deallocate (n_in_cell, p_in_cell)

    end subroutine grid_sweep_3d_other

    !> @brief A helper subroutine to perform a pair search within a cell.
    !> @param cutoff The cutoff distance to find pairs of particles within.
    !> @param ndims The number of spatial dimensions. Used to size input and tmp arrays.
    !> @param xi The LHS particle's coordinates.
    !> @param n The number of RHS particles. Used to size input arrays.
    !> @param x_rhs The RHS particles' coordinates.
    !> @param n_in_cell The number of particles in the given cell.
    !> @param p_in_cell The array of particle indices in the given cell.
    !> @param kernel The SPH kernel to calculate values and gradient values with.
    !> @param pairs The particle_pairs instance to populate with the search.
    !> @param i the LHS particle index.
    pure subroutine sweep_cell(cutoff, ndims, xi, n, x_rhs, n_in_cell, p_in_cell, kernel, pairs, i)

        integer, intent(in):: ndims, n, n_in_cell, p_in_cell(n_in_cell)
        real(fp), intent(in):: cutoff, xi(ndims), x_rhs(ndims, n)
        class(grasph_base_kernel), intent(in):: kernel
        type(particle_pairs), intent(inout):: pairs
        integer, intent(in):: i
        integer:: j, pic
        real(fp):: dx(ndims)

        do pic = 1, n_in_cell
            j = p_in_cell(pic)
            dx(:) = xi(:) - x_rhs(:, j)
            if (sum(dx*dx) < cutoff*cutoff) then
                pairs%npairs_total = pairs%npairs_total + 1
                pairs%pair_ij(1, pairs%npairs_total) = i
                pairs%pair_ij(2, pairs%npairs_total) = j
                call kernel%values(dx, pairs%w(pairs%npairs_total), pairs%dwdx(:, pairs%npairs_total))
            end if
        end do

    end subroutine sweep_cell

end module grasph_pairs
