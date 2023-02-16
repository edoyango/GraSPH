module ORB_sr_m

   use datatypes, only: particles
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
      integer, intent(inout):: ntotal_loc
      type(particles), intent(inout):: parts(:)
      type(neighbour_data), intent(inout):: neighbours(:)
      integer, intent(out):: n_recv_all
      integer:: searchrange(2), entrydepth=0
      real(f):: xmin_loc(dim), xmax_loc(dim)
      logical:: diffuse = .true.
      integer, allocatable:: removal_list(:)

      ! Initialization
      searchrange(:) = [1, ntotal_loc]
      xmin_loc(:) = bounds_loc(1:dim)
      xmin_loc(:) = bounds_loc(dim+1:2*dim)

   end subroutine ORB_sendrecv_diffuse

end module ORB_sr_m