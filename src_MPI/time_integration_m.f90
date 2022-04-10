module time_integration_m

	use globvar,		 only: parts,ntotal_loc,time,cputime,output_time,t_graph,t_dist,test_time,itimestep,maxtimestep,&
        print_step,save_step,maxnloc
	use globvar_para,	 only: procid,numprocs
	use mpi
	use param,			 only: f,dim,dt,tenselem
	
    use input_m,         only: gind,generate_ghost_part
	use flink_list_m,	 only: flink_list
	use ORB_m,			 only: ORB
	use output_m,		 only: output
	use single_step_m,	 only: single_step
    use stress_update_m, only: stress_update
	use summary_m,		 only: print_loadbalance
	
contains
	
	!==============================================================================================================================
	subroutine time_integration
	! Subroutine responsible for the main time-integration loop
	
		implicit none     
		integer:: i,j,k,d,n
		real(f),allocatable:: v_min(:,:),rho_min(:),strain_min(:,:),sig_min(:,:),dvxdt(:,:,:),drho(:,:),dstrain(:,:,:)
        real(f):: dstraini(tenselem)
		
		allocate(v_min(dim,maxnloc),rho_min(maxnloc),strain_min(tenselem,maxnloc),sig_min(tenselem,maxnloc),&
            dvxdt(dim,maxnloc,4),drho(maxnloc,4),dstrain(tenselem,maxnloc,4))
        allocate(gind(maxnloc))
		
		! Time-integration (Leap-Frog)
		do itimestep = 1, maxtimestep
		
			cputime = cputime - MPI_WTIME()
			
			! distributing particles
			call ORB
            
            ! generating ghost particles
            call generate_ghost_part
			
			! Storing velocity and density at initial time-step
			do i = 1,ntotal_loc
				v_min(:,i) = parts(i)%vx(:)
				rho_min(i) = parts(i)%rho
                strain_min(:,i) = parts(i)%strain(:)
                sig_min(:,i) = parts(i)%sig(:)
			end do
			
			!Interaction parameters, calculating neighboring particles
			call flink_list
			
			! calculating forces (k1)
			call single_step(1,dvxdt(:,:,1),drho(:,1),dstrain(:,:,1))
	
			! updating data for mid-timestep base on k1
			do i = 1,ntotal_loc
				parts(i)%vx(:) = v_min(:,i) + 0.5_f*dt*dvxdt(:,i,1)
				parts(i)%rho = rho_min(i) + 0.5_f*dt*drho(i,1)
                dstraini(:) = 0.5_f*dt*dstrain(:,i,1)
                parts(i)%strain(:) = strain_min(:,i) + dstraini(:)
                call stress_update(1,dstraini,sig_min(:,i),parts(i))
			end do
			
			! calculating forces (k2)
			call single_step(2,dvxdt(:,:,2),drho(:,2),dstrain(:,:,2))
	
			! updating data for mid-timestep base on k2
			do i = 1,ntotal_loc
				parts(i)%vx(:) = v_min(:,i) + 0.5_f*dt*dvxdt(:,i,2)
				parts(i)%rho = rho_min(i) + 0.5_f*dt*drho(i,2)
                dstraini(:) = 0.5_f*dt*dstrain(:,i,2)
                parts(i)%strain(:) = strain_min(:,i) + dstraini(:)
                call stress_update(2,dstraini,sig_min(:,i),parts(i))
			end do
			
			! calculating forces (k3)
			call single_step(3,dvxdt(:,:,3),drho(:,3),dstrain(:,:,3))
	
			! updating data for mid-timestep base on k3
			do i = 1,ntotal_loc
				parts(i)%vx(:) = v_min(:,i) + dt*dvxdt(:,i,3)
				parts(i)%rho = rho_min(i) + dt*drho(i,3)
                dstraini(:) = dt*dstrain(:,i,3)
                parts(i)%strain(:) = strain_min(:,i) + dstraini(:)
                call stress_update(3,dstraini,sig_min(:,i),parts(i))
			end do
			
			call single_step(4,dvxdt(:,:,4),drho(:,4),dstrain(:,:,4))
	
			! updating data for mid-timestep base on k1, k2, k3, k4
			do i = 1,ntotal_loc
				parts(i)%vx(:) = v_min(:,i) + dt/6._f*(dvxdt(:,i,1) + 2._f*dvxdt(:,i,2) + 2._f*dvxdt(:,i,3) + dvxdt(:,i,4))
						
				parts(i)%rho = rho_min(i) + dt/6._f*(drho(i,1) + 2._f*drho(i,2) + 2._f*drho(i,3) + drho(i,4))
					
				parts(i)%x(:) = parts(i)%x(:) + dt*parts(i)%vx(:)
                
                dstraini(:) = dt/6._f*(dstrain(:,i,1) + 2._f*dstrain(:,i,2) + 2._f*dstrain(:,i,3) + dstrain(:,i,4))
				parts(i)%strain(:) = strain_min(:,i) + dstraini(:)
                
                call stress_update(4,dstraini,sig_min(:,i),parts(i))
                
			end do
			
			time = time + dt
			
			cputime = cputime + MPI_WTIME()
			
			output_time = output_time - MPI_WTIME()
			
			! write output data
			if (mod(itimestep,save_step).eq.0) then
				call output
			end if

			output_time = output_time + MPI_WTIME()
			
			if (mod(itimestep,print_step).eq.0) then
				call print_loadbalance
			end if
			
		enddo
	
	end subroutine time_integration

end module time_integration_m
