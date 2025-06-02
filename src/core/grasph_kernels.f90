module grasph_kernels

    use grasph_constants, only: fp, pi

    implicit none

    ! private

    type, abstract:: grasph_base_kernel
        integer:: d = 0
        real(fp):: alpha = 0._fp, h = 0._fp, cutoff = 0._fp ! constant smoothing length (for now)
        logical:: initialized = .false.
    contains
        procedure(kernel_init_interface), deferred:: init
        procedure:: values
        procedure(w_interface), deferred:: w, gradw
    end type grasph_base_kernel

    interface
        real(fp) pure function w_interface(self, q)
            import:: fp, grasph_base_kernel
            class(grasph_base_kernel), intent(in):: self
            real(fp), intent(in):: q
        end function w_interface
        subroutine kernel_init_interface(self, d, h)
            import:: fp, grasph_base_kernel
            class(grasph_base_kernel), intent(inout):: self
            integer, intent(in):: d
            real(fp), intent(in):: h
        end subroutine kernel_init_interface
    end interface

    type, extends(grasph_base_kernel):: grasph_cubic_bspline_kernel
    contains
        procedure:: init => init_cubic_bspline, w => w_cubic_bspline, gradw => gradw_cubic_bspline
    end type grasph_cubic_bspline_kernel

    public:: grasph_base_kernel, grasph_cubic_bspline_kernel
contains
    pure subroutine values(self, dx, w, dwdx)
        class(grasph_base_kernel), intent(in):: self
        real(fp), intent(in):: dx(self%d)
        real(fp), intent(out):: w, dwdx(self%d)
        real(fp):: q, r
        r = sqrt(sum(dx**2))
        q = r/self%h
        w = self%w(q)
        dwdx(:) = self%gradw(q)*dx(:)/(r*self%h)
    end subroutine values

    subroutine init_cubic_bspline(self, d, h)
        class(grasph_cubic_bspline_kernel), intent(inout):: self
        integer, intent(in):: d
        real(fp), intent(in):: h
        self%d = d
        self%h = h
        self%cutoff = 2._fp*h
        if (d == 2) then
            self%alpha = 10._fp/(7._fp*pi*h*h)
        elseif (d == 3) then
            self%alpha = 1._fp/(pi*h*h*h)
        endif
        self%initialized = .true.
    end subroutine init_cubic_bspline

    real(fp) pure function w_cubic_bspline(self, q)
        class(grasph_cubic_bspline_kernel), intent(in):: self
        real(fp), intent(in):: q
        w_cubic_bspline = self%alpha*(0.25_fp*dim(2._fp, q)**3 - dim(1._fp, q)**3)
    end function w_cubic_bspline

    real(fp) pure function gradw_cubic_bspline(self, q)
        class(grasph_cubic_bspline_kernel), intent(in):: self
        real(fp), intent(in):: q
        gradw_cubic_bspline = -self%alpha*3._fp*(0.25_fp*dim(2._fp, q)**2 - dim(1._fp, q)**2)
    end function gradw_cubic_bspline

end module grasph_kernels