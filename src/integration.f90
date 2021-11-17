!=======================================================================================================
subroutine time_integration( )
!=======================================================================================================

use param
implicit none     
integer:: i,j,k,d,n
real(8):: t1,t2
real(8),allocatable:: v_min(:,:),rho_min(:),dvxdt(:,:,:),drho(:,:)
allocate(v_min(dim,ntotal),rho_min(ntotal),dvxdt(dim,ntotal,4),drho(ntotal+nvirt,4))
testtime(:) = 0d0
call CPU_TIME(s1)

! Time-integration (Leap-Frog)
do itimestep = nstart+1, nstart+maxtimestep   
	   
  current_ts=current_ts+1
  if (mod(itimestep,print_step).eq.0) then
	write(*,*)'______________________________________________'
	write(*,*)'  current number of time step =', itimestep,'     current time=', real(time+dt)
	write(*,*)'______________________________________________'
  endif      
  
  ! Update density and velocity half a time step (not at first time-step)
  dvxdt(:,:,:) = 0d0
  drho(:,:) = 0d0
  do i = 1,ntotal
    v_min(:,i) = parts(i)%vx(:)
    rho_min(i) = parts(i)%rho
  end do
  
  call CPU_TIME(t1)
  !Interaction parameters, calculating neighboring particles
  call flink_list
  call CPU_TIME(t2)
  testtime(1) = testtime(1) + t2 - t1
  ! calculating forces (k1)
  call single_step(1,dvxdt(1:dim,1:ntotal,1),drho(1:ntotal,1))
  
  ! updating data for mid-timestep base on k1
  do i = 1,ntotal
    parts(i)%vx(:) = v_min(:,i) + 0.5d0*dt*dvxdt(:,i,1)
    parts(i)%rho = rho_min(i) + 0.5d0*dt*drho(i,1)
  end do
  
  ! calculating forces (k2)
  call single_step(2,dvxdt(1:dim,1:ntotal,2),drho(1:ntotal,2))
  
  ! updating data for mid-timestep base on k2
  do i = 1,ntotal
    parts(i)%vx(:) = v_min(:,i) + 0.5d0*dt*dvxdt(:,i,2)
    parts(i)%rho = rho_min(i) + 0.5d0*dt*drho(i,2)
  end do
  
  ! calculating forces (k3)
  call single_step(3,dvxdt(1:dim,1:ntotal,3),drho(1:ntotal,3))
  
  ! updating data for mid-timestep base on k3
  do i = 1,ntotal
    parts(i)%vx(:) = v_min(:,i) + dt*dvxdt(:,i,3)
    parts(i)%rho = rho_min(i) + dt*drho(i,3)
  end do
  
  call single_step(4,dvxdt(1:dim,1:ntotal,4),drho(1:ntotal,4))
  
  ! updating data for mid-timestep base on k1, k2, k3, k4
  do i = 1,ntotal
    parts(i)%vx(:) = v_min(:,i) + dt/6d0*(dvxdt(:,i,1) + 2d0*dvxdt(:,i,2) + 2d0*dvxdt(:,i,3) + dvxdt(:,i,4))
	                 
    parts(i)%rho = rho_min(i) + dt/6d0*(drho(i,1) + 2d0*drho(i,2) + 2d0*drho(i,3) + drho(i,4))
	               
	parts(i)%x(:) = parts(i)%x(:) + dt*parts(i)%vx(:)
  end do

  time = time + dt
  call CPU_TIME(s2)

  ! write output data
  if (mod(itimestep,save_step).eq.0) then
	n = itimestep / save_step
    call output( )
  endif 

enddo

nstart=current_ts
cputime = cputime+s2-s1
write(*,*) testtime(1:5)
end
