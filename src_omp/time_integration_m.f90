module time_integration_m

	use globvar, only: time,cputime,output_time,test_time,ntotal,parts,print_step,save_step,itimestep,maxtimestep
	use param, only: f,dim,rh0,gamma,c,dt
	use omp_lib
	
	use flink_list_m, only: flink_list
	use output_m, only: output
	use single_step_m, only: single_step
	use summary_m, only: print_update
	
	public:: time_integration

contains

	!==============================================================================================================================
	subroutine time_integration( )
	! Subroutine responsible for the main time-integration loop
	
		implicit none     
		integer:: i,j,k,d,n
		real(f):: t1,t2
		real(f),allocatable:: v_min(:,:),rho_min(:),dvxdt(:,:,:),drho(:,:)
		
		allocate(v_min(dim,ntotal),rho_min(ntotal),dvxdt(dim,ntotal,4),drho(ntotal,4))
		
		! Time-integration (Leap-Frog)
		do itimestep = 1, maxtimestep
		
			t1 = OMP_GET_WTIME()
			
			!Interaction parameters, calculating neighboring particles
			call flink_list
			
			! Update density and velocity half a time step (not at first time-step)
			!$OMP PARALLEL DO DEFAULT(SHARED)
			do i = 1,ntotal
				v_min(:,i) = parts(i)%vx(:)
				rho_min(i) = parts(i)%rho
				parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma-1d0)/gamma
			end do
			!$OMP END PARALLEL DO
			
			! calculating forces (k1)
			call single_step(1,dvxdt(:,:,1),drho(:,1))
			
			! updating data for mid-timestep base on k1
			!$OMP PARALLEL DO DEFAULT(SHARED)
			do i = 1,ntotal
				parts(i)%vx(:) = v_min(:,i) + 0.5d0*dt*dvxdt(:,i,1)
				parts(i)%rho = rho_min(i) + 0.5d0*dt*drho(i,1)
				parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma-1d0)/gamma
			end do
			!$OMP END PARALLEL DO
			
			! calculating forces (k2)
			call single_step(2,dvxdt(:,:,2),drho(:,2))
			
			! updating data for mid-timestep base on k2
			!$OMP PARALLEL DO DEFAULT(SHARED)
			do i = 1,ntotal
				parts(i)%vx(:) = v_min(:,i) + 0.5d0*dt*dvxdt(:,i,2)
				parts(i)%rho = rho_min(i) + 0.5d0*dt*drho(i,2)
				parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma-1d0)/gamma
			end do
			!$OMP END PARALLEL DO
			
			! calculating forces (k3)
			call single_step(3,dvxdt(:,:,3),drho(:,3))
			
			! updating data for mid-timestep base on k3
			!$OMP PARALLEL DO DEFAULT(SHARED)
			do i = 1,ntotal
				parts(i)%vx(:) = v_min(:,i) + dt*dvxdt(:,i,3)
				parts(i)%rho = rho_min(i) + dt*drho(i,3)
				parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma-1d0)/gamma
			end do
			!$OMP END PARALLEL DO
			
			call single_step(4,dvxdt(:,:,4),drho(:,4))
			
			! updating data for mid-timestep base on k1, k2, k3, k4
			!$OMP PARALLEL DO DEFAULT(SHARED)
			do i = 1,ntotal
				parts(i)%vx(:) = v_min(:,i) + dt/6d0*(dvxdt(:,i,1) + 2d0*dvxdt(:,i,2) + 2d0*dvxdt(:,i,3) + dvxdt(:,i,4))
								
				parts(i)%rho = rho_min(i) + dt/6d0*(drho(i,1) + 2d0*drho(i,2) + 2d0*drho(i,3) + drho(i,4))
							
				parts(i)%x(:) = parts(i)%x(:) + dt*parts(i)%vx(:)
				
				parts(i)%p = rh0*c**2*((parts(i)%rho/rh0)**gamma-1d0)/gamma
			end do
			!$OMP END PARALLEL DO
			
			time = time + dt
			
			t2 = OMP_GET_WTIME()
			cputime = cputime + t2 - t1
			
			t1 = OMP_GET_WTIME()
			
			! write output data
			if (mod(itimestep,save_step).eq.0) then
				call output( )	
			end if 
			
			t2 = OMP_GET_WTIME()
			output_time = output_time + t2 - t1
	
			if (mod(itimestep,print_step).eq.0) then
				call print_update
			end if
		
		enddo
		
	end subroutine time_integration

end module time_integration_m