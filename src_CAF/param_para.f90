! paramaters that relate to the parallel scheme for this code e.g. array sizes, subdomain boundary
! update frequency, number of msgs etc.
module param_para

   use datatypes, only: particles
   use param, only: dim, f

   implicit none
   private
   ! parameters for ORB algorithm
! dcell_ORB = size of each cell used to draw subdomain boundaries (in multiples of hsml)
!             smaller -> better partition, but slower and requires more memory.
!             bigger -> poorer partition, but faster and requires less memory.
! bound_extend = the distance of which to extend the current global domain boundaries (in multiples of kernel radii). Choose a
!                large value if you expect global extents of your simulation to change a lot over the simulation.
! box_ratio_threshold = controls the required proportional change in aspect ratio before subdomain bounadry orientations are
!                       re-determined.
!                       Smaller value -> more frequent axes reorientations of cuts.
! ORBcheck1, ORBcheck2 = how frequent to check if boundaries need updating. Lets say the boundar is updated at time-step N, the
!                        next check will be at time-step N + ORBcheck1. If the boundaries are not updated, checks will occur every
!                        ORBcheck2 time-steps thereafter. Too small and program will be very slow. Too large, and load balance will
!                        be very poor in general. The more slowly the particles move, ORBchecks can be less frequent.
   real(f), parameter, public:: dcell_ORB = 1_f, box_ratio_threshold = 0.25_f, bound_extend = 10_f
   integer, parameter, public:: ORBcheck1 = 200, ORBcheck2 = 100

   ! type to hold variables related to partitioning
   type, public:: partition_tracking
      integer:: mintstep_bn_part = HUGE(1), mintstep_bn_reorient = HUGE(1), maxtstep_bn_part = 0, &
                maxtstep_bn_reorient = 0, prev_part_tstep, prev_reorient_tstep, n_parts = 0, n_reorients = 0
   end type partition_tracking

   ! type to hold all variables related to data transfer between neighbouring processes
   type, public:: neighbour_data
      integer:: pid, nphys_send, nphys_recv, nhalo_send, nhalo_recv
      real(f):: bounds(2*dim)
      integer, allocatable:: halo_pindex(:)
   end type neighbour_data

end module param_para
