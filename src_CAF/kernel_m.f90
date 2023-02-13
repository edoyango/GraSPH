module kernel_m

   use param, only: f

   private
   public:: kernel_k

contains

   !====================================================================================================================
   pure function kernel_k(skf) result(scale_k)
   ! setting k parameter for kernel radius (kh)

      implicit none
      integer, intent(in):: skf
      real(f):: scale_k

      select case(skf)
      case default
         scale_k = 2._f
      case(2)
         scale_k = 2.5_f
      case(3, 7)
        scale_k = 3._f
      end select

    end function kernel_k

end module kernel_m