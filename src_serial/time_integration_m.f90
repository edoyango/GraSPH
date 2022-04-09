module time_integration_m

	use globvar, only: time,cputime,output_time,test_time,ntotal,nvirt,parts,print_step,save_step,itimestep,maxtimestep
    use input_m, only: gind
	use param, only: f,dim,dt,tenselem
	
	use flink_list_m, only: flink_list
    use input_m, only: generate_ghost_part,update_ghost_part,update_virt_part
	use output_m, only: output
	use single_step_m, only: single_step
    use stress_update_m, only: stress_update
	use summary_m, only: print_update
	
	public:: time_integration

contains

	!==============================================================================================================================
	subroutine time_integration( )
	! Subroutine responsible for the main time-integration loop
	
		implicit none     
		integer:: i,j,k,d,n
		real(f):: t1,t2,t3,t4,dstraini(tenselem)
		real(f),allocatable:: v_min(:,:),rho_min(:),strain_min(:,:),sig_min(:,:),dvxdt(:,:,:),drho(:,:),dstrain(:,:,:)
		
		allocate(v_min(dim,ntotal),rho_min(ntotal),dvxdt(dim,ntotal,4),drho(ntotal,4))
        allocate(strain_min(tenselem,ntotal),dstrain(tenselem,ntotal,4),sig_min(tenselem,ntotal))
        allocate(gind(ntotal))
        
        call CPU_TIME(t1)
		
		! Time-integration (Leap-Frog)
		do itimestep = 1, maxtimestep
            
            ! Update density and velocity half a time step (not at first time-step)
			do i = 1,ntotal
				v_min(:,i) = parts(i)%vx(:)
				rho_min(i) = parts(i)%rho
                strain_min(:,i) = parts(i)%strain(:)
                sig_min(:,i) = parts(i)%sig(:)
			end do
            
            call generate_ghost_part
			
			!Interaction parameters, calculating neighboring particles
			call flink_list
            
			! calculating forces (k1)
			call single_step(1,dvxdt(:,:,1),drho(:,1),dstrain(:,:,1))
			
			! updating data for mid-timestep base on k1
			do i = 1,ntotal
				parts(i)%vx(:) = v_min(:,i) + 0.5_f*dt*dvxdt(:,i,1)
				parts(i)%rho = rho_min(i) + 0.5_f*dt*drho(i,1)
                dstraini = 0.5_f*dt*dstrain(:,i,1)
                parts(i)%strain(:) = strain_min(:,i) + dstraini(:)
                call stress_update(1,dstraini,sig_min(:,i),parts(i))
			end do
            
            call update_ghost_part
            
			! calculating forces (k2)
			call single_step(2,dvxdt(:,:,2),drho(:,2),dstrain(:,:,2))
			
			! updating data for mid-timestep base on k2
			do i = 1,ntotal
				parts(i)%vx(:) = v_min(:,i) + 0.5_f*dt*dvxdt(:,i,2)
				parts(i)%rho = rho_min(i) + 0.5_f*dt*drho(i,2)
                dstraini = 0.5_f*dt*dstrain(:,i,2)
                parts(i)%strain(:) = strain_min(:,i) + dstraini(:)
                call stress_update(2,dstraini,sig_min(:,i),parts(i))
			end do
            
            call update_ghost_part
            
			! calculating forces (k3)
			call single_step(3,dvxdt(:,:,3),drho(:,3),dstrain(:,:,3))
			
			! updating data for mid-timestep base on k3
			do i = 1,ntotal
				parts(i)%vx(:) = v_min(:,i) + dt*dvxdt(:,i,3)
				parts(i)%rho = rho_min(i) + dt*drho(i,3)
                dstraini = dt*dstrain(:,i,3)
                parts(i)%strain(:) = strain_min(:,i) + dstraini
                call stress_update(3,dstraini,sig_min(:,i),parts(i))
			end do
            
            call update_ghost_part
            
			call single_step(4,dvxdt(:,:,4),drho(:,4),dstrain(:,:,4))
			
			! updating data for mid-timestep base on k1, k2, k3, k4
			do i = 1,ntotal
				parts(i)%vx(:) = v_min(:,i) + dt/6._f*(dvxdt(:,i,1) + 2._f*dvxdt(:,i,2) + 2._f*dvxdt(:,i,3) + dvxdt(:,i,4))
								
				parts(i)%rho = rho_min(i) + dt/6._f*(drho(i,1) + 2._f*drho(i,2) + 2._f*drho(i,3) + drho(i,4))
							
				parts(i)%x(:) = parts(i)%x(:) + dt*parts(i)%vx(:)
                
                dstraini = dt/6._f*(dstrain(:,i,1) + 2._f*dstrain(:,i,2) + 2._f*dstrain(:,i,3) + dstrain(:,i,4))
                parts(i)%strain(:) = strain_min(:,i) + dstraini(:)
				
                call stress_update(4,dstraini,sig_min(:,i),parts(i))
			end do
            
            call update_ghost_part
            
			time = time + dt
			
			call CPU_TIME(t3)
			
			! write output data
			if (mod(itimestep,save_step).eq.0) then
				call output( )
			end if 
			
			call CPU_TIME(t4)
			output_time = output_time + t4 - t3
	
			if (mod(itimestep,print_step).eq.0) then
                call CPU_TIME(t2)
                cputime = t2 - t1
				call print_update
			end if
		
		enddo
    
        call CPU_TIME(t2)
        cputime = t2 - t1
		
	end subroutine time_integration

end module time_integration_m
