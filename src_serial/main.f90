program SPH

   use datatypes, only: particles, interactions
   use param, only: skf, f

   use input_m, only: generate_real_part, generate_virt_part, allocatePersistentArrays, generate_ghost_part, return_ntotal, &
                      return_nvirt
   use kernel_m, only: kernel_k
   use output_m, only: output
   use summary_m, only: preamble, time_print, print_summary
   use time_integration_m, only: time_integration

   implicit none
   type(particles), allocatable:: parts(:) ! particle array
   type(interactions), allocatable:: pairs(:) ! interaction array
   integer:: ntotal, nvirt, nghos, niac ! tracking no. of particles, and interactions
   integer:: maxtimestep, save_step, print_step ! timestep related variables
   real(f):: time = 0_f, scale_k
   real(f):: cputime = 0_f, output_time = 0_f, test_time = 0_f !measuring compute time
   integer, allocatable:: gind(:)

   !Printing preamble to screen
   call preamble(maxtimestep, print_step, save_step)

   ! setting k parameter for kernel radius (r = kh)
   scale_k = kernel_k(skf)

   ntotal = return_ntotal()
   nvirt = return_nvirt()

   write (*, '(A24,1x,I9,1x,A19)') 'Total simulation size of', ntotal, 'physical particles.'
   write (*, '(A24,1x,I9,1x,A19)') '                        ', nvirt, 'virtual particles.'

   call allocatePersistentArrays(ntotal, nvirt, parts, pairs, gind)

   !Creat physical and virtual boundary particles
   call generate_real_part(ntotal, parts)
   call generate_virt_part(ntotal, nvirt, parts)
   call generate_ghost_part(scale_k, ntotal, nvirt, nghos, parts, gind)
   call output(0, save_step, ntotal, nvirt, nghos, parts)

   !Entering discretized time-integration loop
   call time_integration(time, cputime, output_time, test_time, ntotal, nvirt, nghos, parts, print_step, save_step, &
                         maxtimestep, niac, pairs, scale_k, gind)

   !Printing post-amble to terminal
   call time_print
   call print_summary(cputime, output_time)

end
