module grasph_pairs

    use grasph_constants, only: fp
    use grasph_kernels, only: grasph_base_kernel

    implicit none

    private

    type particle_pairs
        integer:: n = 0, npairs_per_particle = 0, npairs_total = 0
        integer, allocatable:: rhs(:), offsets(:)
        real(fp), allocatable:: w(:), dwdx(:, :)
        logical:: initialized = .false.
    contains
        procedure:: init => particle_pairs_init
    end type particle_pairs

    ! abstract interface provided to be used in function pointers
    abstract interface
        type(particle_pairs) pure function find_pairs_fixed_h(x, ndims, n, cutoff, kernel, npairs_per_particle)
            import:: fp, grasph_base_kernel, particle_pairs
            integer, intent(in):: ndims, n, npairs_per_particle
            real(fp), intent(in):: x(ndims, n), cutoff
            class(grasph_base_kernel), intent(in):: kernel
        end function find_pairs_fixed_h
    end interface

    public:: particle_pairs, dsearch, cell_list_search, find_pairs_fixed_h

contains

    pure subroutine particle_pairs_init(self, n, npairs_per_particle, ndims)

        class(particle_pairs), intent(inout):: self
        integer, intent(in):: n, npairs_per_particle, ndims

        self%n = n
        self%npairs_per_particle = npairs_per_particle
        self%npairs_total = 0

        if (self%initialized) deallocate(self%rhs, self%offsets, self%w, self%dwdx)
        allocate(self%rhs(n*npairs_per_particle))
        allocate(self%offsets(n+1), source=0)
        allocate(self%w(n*npairs_per_particle), self%dwdx(ndims, n*npairs_per_particle))

        self%initialized = .true.

    end subroutine particle_pairs_init

    type(particle_pairs) pure function dsearch(x, ndims, n, cutoff, kernel, npairs_per_particle) result(pairs)

        integer, intent(in):: ndims, n, npairs_per_particle
        real(fp), intent(in):: x(ndims, n), cutoff
        class(grasph_base_kernel), intent(in):: kernel
        integer:: i, j
        real(fp):: dx(ndims)

        call pairs%init(n, npairs_per_particle, ndims)

        do i = 1, n-1
            do j = i+1, n
                dx(:) = x(:, i) - x(:, j)
                if (sum(dx(:)**2) < cutoff*cutoff) then
                    pairs%npairs_total = pairs%npairs_total + 1
                    pairs%rhs(pairs%npairs_total) = j
                    call kernel%values(dx, pairs%w(pairs%npairs_total), pairs%dwdx(:, pairs%npairs_total))
                endif
            enddo
            pairs%offsets(i+1) = pairs%npairs_total
        enddo
        pairs%offsets(n+1) = pairs%npairs_total
    end function dsearch

    type(particle_pairs) pure function cell_list_search(x, ndims, n, cutoff, kernel, npairs_per_particle) result(pairs)

        integer, intent(in):: ndims, n, npairs_per_particle
        real(fp), intent(in):: x(ndims, n), cutoff
        class(grasph_base_kernel), intent(in):: kernel
        real(fp):: minextents(ndims), maxextents(ndims), dcell
        integer:: i, ngridx(ndims), grid_idx(ndims, n) ! might need to be allocatable in the future...

        ! define grid
        dcell = cutoff ! not functionally meaningful, but helpful conceptually
        minextents(:) = minval(x, 2) - 2._fp*dcell
        maxextents(:) = maxval(x, 2) + 2._fp*dcell
        ngridx(:) = int((maxextents(:)-minextents(:))/dcell) + 1
        ! technically, maxextents should be adjusted, but it isn't used from herin

        do i = 1, n
            grid_idx(:, i) = int((x(:, i) - minextents(:))/dcell) + 1
        enddo

        call pairs%init(n, npairs_per_particle, ndims)

        select case (ndims)
        case(2)
            call grid_sweep_2d(n, cutoff, npairs_per_particle, kernel, ngridx, grid_idx, x, pairs)
        case(3)
            call grid_sweep_3d(n, cutoff, npairs_per_particle, kernel, ngridx, grid_idx, x, pairs)
        case default
            error stop "cell_list_search: only 2d and 3d cases are supported!"
        end select
    
    end function cell_list_search

    pure subroutine grid_sweep_2d(n, cutoff, npairs_per_particle, kernel, ngridx, grid_idx, x, pairs)

        integer, intent(in):: n, npairs_per_particle
        real(fp), intent(in):: cutoff, x(2, n)
        class(grasph_base_kernel), intent(in):: kernel
        integer, intent(in):: ngridx(2), grid_idx(2, n)
        type(particle_pairs), intent(inout):: pairs
        integer:: i, j, icell, jcell, jj, pic
        real(fp):: dx(2)
        integer, allocatable:: n_in_cell(:, :), p_in_cell(:, :, :)

        allocate(n_in_cell(ngridx(1), ngridx(2)), source=0)
        allocate(p_in_cell(npairs_per_particle, ngridx(1), ngridx(2)))

        ! populate grid
        do i = 1, n
            icell = grid_idx(1, i)
            jcell = grid_idx(2, i)
            n_in_cell(icell, jcell) = n_in_cell(icell, jcell) + 1
            p_in_cell(n_in_cell(icell, jcell), icell, jcell) = i
        enddo

        ! sweep
        do i = 1, n
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
                        pairs%rhs(pairs%npairs_total) = j
                        call kernel%values(dx, pairs%w(pairs%npairs_total), pairs%dwdx(:, pairs%npairs_total))
                    endif
                endif
            enddo
            ! right cell
            icell = icell + 1
            call sweep_cell(i, cutoff, 2, n, x, n_in_cell(icell, jcell), &
                            p_in_cell(:, icell, jcell), kernel, pairs)
            ! top row
            jcell = jcell + 1
            do icell = grid_idx(1, i) - 1, grid_idx(1, i) + 1
                call sweep_cell(i, cutoff, 2, n, x, n_in_cell(icell, jcell), &
                                p_in_cell(:, icell, jcell), kernel, pairs)
            enddo
            pairs%offsets(i+1) = pairs%npairs_total
        enddo

        deallocate(n_in_cell, p_in_cell)

    end subroutine grid_sweep_2d

    pure subroutine grid_sweep_3d(n, cutoff, npairs_per_particle, kernel, ngridx, grid_idx, x, pairs)

        integer, intent(in):: n, npairs_per_particle
        real(fp), intent(in):: cutoff, x(3, n)
        class(grasph_base_kernel), intent(in):: kernel
        integer, intent(in):: ngridx(3), grid_idx(3, n)
        type(particle_pairs), intent(inout):: pairs
        integer:: i, j, icell, jcell, kcell, jj, pic
        real(fp):: dx(3)
        integer, allocatable:: n_in_cell(:, :, :), p_in_cell(:, :, :, :)

        allocate(n_in_cell(ngridx(1), ngridx(2), ngridx(3)), source=0)
        allocate(p_in_cell(npairs_per_particle, ngridx(1), ngridx(2), ngridx(3)))

        ! populate grid
        n_in_cell(:, :, :) = 0
        do i = 1, n
            icell = grid_idx(1, i)
            jcell = grid_idx(2, i)
            kcell = grid_idx(3, i)
            n_in_cell(icell, jcell, kcell) = n_in_cell(icell, jcell, kcell) + 1
            p_in_cell(n_in_cell(icell, jcell, kcell), icell, jcell, kcell) = i
        enddo

        ! sweep
        do i = 1, n
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
                        pairs%rhs(pairs%npairs_total) = j
                        call kernel%values(dx, pairs%w(pairs%npairs_total), pairs%dwdx(:, pairs%npairs_total))
                    endif
                endif
            enddo
            ! right cell
            icell = icell + 1
            call sweep_cell(i, cutoff, 3, n, x, n_in_cell(icell, jcell, kcell), &
                            p_in_cell(:, icell, jcell, kcell), kernel, pairs)
            ! north-middle layer
            jcell = jcell + 1
            do icell = grid_idx(1, i) - 1, grid_idx(1, i) + 1
                call sweep_cell(i, cutoff, 3, n, x, n_in_cell(icell, jcell, kcell), &
                                p_in_cell(:, icell, jcell, kcell), kernel, pairs)
            enddo
            ! top layer
            kcell = kcell + 1
            do jcell = grid_idx(2, i) - 1, grid_idx(2, i) + 1
                do icell = grid_idx(1, i) - 1, grid_idx(1, i) + 1
                    call sweep_cell(i, cutoff, 3, n, x, n_in_cell(icell, jcell, kcell), &
                                    p_in_cell(:, icell, jcell, kcell), kernel, pairs)
                enddo
            enddo
            pairs%offsets(i+1) = pairs%npairs_total
        enddo

        deallocate(n_in_cell, p_in_cell)

    end subroutine grid_sweep_3d

    pure subroutine sweep_cell(i, cutoff, ndims, n, x, n_in_cell, p_in_cell, kernel, pairs)

        integer, intent(in):: i, ndims, n, n_in_cell, p_in_cell(n_in_cell)
        real(fp), intent(in):: cutoff, x(ndims, n)
        class(grasph_base_kernel), intent(in):: kernel
        type(particle_pairs), intent(inout):: pairs
        integer:: j, pic
        real(fp):: dx(ndims)

        do pic = 1, n_in_cell
            j = p_in_cell(pic)
            dx(:) = x(:, i) - x(:, j)
            if (sum(dx*dx) < cutoff*cutoff) then
                pairs%npairs_total = pairs%npairs_total + 1
                pairs%rhs(pairs%npairs_total) = j
                call kernel%values(dx, pairs%w(pairs%npairs_total), pairs%dwdx(:, pairs%npairs_total))
            endif
        enddo

    end subroutine sweep_cell

end module grasph_pairs