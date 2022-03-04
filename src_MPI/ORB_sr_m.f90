module ORB_sr_m

	use globvar, 		only: ntotal_loc,nhalo_loc,parts,maxnloc,itimestep,scale_k
	use globvar_para,	only: procid,numprocs,PhysPackSend,bounds_glob,n_process_neighbour,proc_neighbour_list,nphys_send,&
						nphys_recv,parttype,ierr
	use mpi
	use param,			only: f,dim,hsml
	use error_msg_m,	only: error_msg
	
	public:: ORB_sendrecv_diffuse,ORB_sendrecv_halo,ORB_sendrecv_haloupdate
	
contains

	!==============================================================================================================================
	recursive function ORB_sendrecv_diffuse( entrydepth,searchrange,nrequest,request,n_recv_all ) result (success)
	! Recursive function to exchange physical particles. In cases were subdomain boundaries are updated, the possibility of needing
	! diffusion is considered
	
		use globvar_para,	only: repartition_mode
		
		implicit none
		integer,intent(in):: searchrange(2)
		integer,intent(inout):: nrequest,request(2*n_process_neighbour),n_recv_all,entrydepth
		integer:: d,i,j,k,pid,n,pos_recv,success
		integer:: removal_list(searchrange(2)-searchrange(1)+2),nphys_send_all,diff_dest,ndiffuse_loc,ndiffuse_all,searchrange_next(2)
		real(f):: xmin_loc(dim),xmax_loc(dim),xmin_rem(dim),xmax_rem(dim),xi(dim),dr,dr_min
		integer:: status(2*n_process_neighbour,MPI_STATUS_SIZE)
		
		! Initialization
		success = 1
		
		! if entrydepth > 1, then diffusion is occuring. Wait for ongoing comms to complete. Update ntotal_loc as required.
		if (entrydepth.gt.0) then
			call MPI_WAITALL(nrequest,request,status,ierr)
			ntotal_loc = ntotal_loc + n_recv_all
			deallocate( PhysPackSend )
		end if
		
		xmin_loc(:) = bounds_glob(1:dim,procid+1)
		xmax_loc(:) = bounds_glob(dim+1:2*dim,procid+1)
			
		! Searching particles to remove within indices of searchrange(1) and searchrange(2), inclusive.
		! At node 0, searchrange(1:2) = [1,ntotal_loc)
		nphys_send_all = 0
		do i = searchrange(1),searchrange(2)
			xi(:) = parts(i)%x(:)
			if ( any( [xi(:).lt.xmin_loc(:) , xi(:).ge.xmax_loc(:)] ) ) then
				nphys_send_all = nphys_send_all + 1
				removal_list(nphys_send_all) = i
			end if
		end do
		removal_list(nphys_send_all+1) = 0 ! This is needed due to a quirk in the loop below
		
		! If there are any particles that do not belong to the host process,
		! begin searching for neighbouring processes to send the particle to.
		! If particle is not contained with subdomain boundaries, send the particle to nearest process neighbouring current host.
		allocate( PhysPackSend(nphys_send_all,n_process_neighbour) )
		nphys_send(:) = 0
		ndiffuse_loc = 0
		if (nphys_send_all.gt.0) then
			loop_through_parts: do j = 1,nphys_send_all
				i = removal_list(j)
				xi(:) = parts(i)%x(:) !placeholder variable of particle position for code readibility
				do n = 1,n_process_neighbour
					pid = proc_neighbour_list(n) + 1
					xmin_rem(:) = bounds_glob(1:dim,pid) ! min boundary for remote process, pid
					xmax_rem(:) = bounds_glob(dim+1:2*dim,pid) ! max boundary for remote process, pid
					if ( all( [xi(:).ge.xmin_rem(:) , xi(:).lt.xmax_rem(:)] ) ) then
		
						nphys_send(n) = nphys_send(n) + 1
						PhysPackSend(nphys_send(n),n) = parts(i)
		
						cycle loop_through_parts ! stop searching for a subdomain neighbour and move to next particle
					end if
				end do
	
				! particle belongs to non-neighbouring process. 
				! Evaluates closeness by considering minimum distance between particle and a processes edge,face, or vertice
				ndiffuse_loc = ndiffuse_loc + 1
				dr_min = huge(1_f)
				do n = 1,n_process_neighbour
					pid = proc_neighbour_list(n) + 1
					dr = 0_f
					do d = 1,dim
						dr = dr + MAX(0d0,bounds_glob(d,pid)-parts(i)%x(d),parts(i)%x(d)-bounds_glob(dim+d,pid))**2
					end do
					dr = SQRT(dr)
					if (dr.lt.dr_min) then
						diff_dest = n
						dr_min = dr
					end if
				end do
	
				nphys_send(diff_dest) = nphys_send(diff_dest) + 1
				PhysPackSend(nphys_send(diff_dest),diff_dest) = parts(i)
				
				!call error_msg(itimestep,procid,4,ind(i))
			end do loop_through_parts
		end if
		
		! Posting non-blocking send/recv to exchange info with neighers of no. particles being sent
		do n = 1,n_process_neighbour
			pid = proc_neighbour_list(n)
			call MPI_IRECV(nphys_recv(n),1,MPI_INTEGER,pid,0,MPI_COMM_WORLD,request(2*n-1),ierr)
			call MPI_ISEND(nphys_send(n),1,MPI_INTEGER,pid,0,MPI_COMM_WORLD,request(2*n),ierr)
		end do
		
		! Shifting information to remove sent particles
		n = 0
		if ( (nphys_send_all.gt.0) .and. (nphys_send_all.lt.ntotal_loc) ) then
			do i = removal_list(1),ntotal_loc
				if ( i.eq.removal_list(n+1) ) then
					n = n + 1
				else
					parts(i-n) = parts(i) ! moving particle to new local position
				end if
			end do
		end if
		
		ntotal_loc = ntotal_loc - nphys_send_all
		
		! Wait for non-blocking send to complete
		call MPI_WAITALL(2*n_process_neighbour,request,status,ierr)
		
		! Calculating total number of physical particles to be received. Check if this exceeds particle array boundaries
		n_recv_all = SUM(nphys_recv)
		if (ntotal_loc+n_recv_all.gt.maxnloc) call error_msg(2,parts(ntotal_loc)%indloc)
		
		! Non-blocking sends to exchange physical particles that have moved processes ---------------------------------------------------
		pos_recv = ntotal_loc + 1
		nrequest = 0
		do n = 1,n_process_neighbour
	
			pid = proc_neighbour_list(n)
	
			if (nphys_recv(n).gt.0) then
	
				nrequest = nrequest + 1
				call MPI_IRECV(parts(pos_recv)%indglob,nphys_recv(n),parttype,pid,0,MPI_COMM_WORLD,request(nrequest),ierr)
				
				pos_recv = pos_recv + nphys_recv(n)
					
			end if
	
			if (nphys_send(n).gt.0) then
			
				nrequest = nrequest + 1
				call MPI_ISEND(PhysPackSend(1,n)%indglob,nphys_send(n),parttype,pid,0,MPI_COMM_WORLD,request(nrequest),ierr)
		
			end if
	
		end do
		
		! if subdomain boundary update has occurred, check if diffusion is necessary.
		! Perform if necessary
		if (repartition_mode.ge.2) then
		
			! Checking if any process has particles that need to be diffused
			call MPI_ALLREDUCE(ndiffuse_loc,ndiffuse_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)					
			
			! If any particles need to be diffused, repeat process on newly received particles
			if (ndiffuse_all.ne.0) then
				
				if (procid.eq.0) then
					write(*,'(A48,I7)') 'Diffusion occuring... Current timestep:         ',itimestep
					write(*,'(A48,I7)') '                      Current depth:            ',entrydepth
					write(*,'(A48,I7)') '                      Particles to be diffused: ',ndiffuse_all
				end if
				
				entrydepth = entrydepth + 1
				searchrange_next = (/ntotal_loc+1,ntotal_loc+n_recv_all/)
				success = ORB_sendrecv_diffuse( entrydepth,searchrange_next,nrequest,request,n_recv_all )
			end if
		end if
		success = 0
		
	end function ORB_sendrecv_diffuse
	
	!==============================================================================================================================
	subroutine ORB_sendrecv_halo( request_in,request_out,nphys_recv_all,nrequest )
	
		!subroutine responsible for sending sending halo particle information between processes, given 
		!predetermined subdomain boundaires. 
		!Note: subdomain boundaries are used as inputs (bounds_glob). 
		
		use globvar_para,	only: halo_pindex,nhalo_send,nhalo_recv,halotype_indexed,haloupdatetype_indexed
		use param_para,		only: halotype,haloupdatetype
	
		implicit none
		integer,intent(inout):: request_in(2*n_process_neighbour)
		integer,intent(inout):: nphys_recv_all,nrequest
		integer,intent(out):: request_out(2*n_process_neighbour)
		integer:: status(MPI_STATUS_SIZE,2*n_process_neighbour)
		integer:: d,i,j,k,pid,n_send_all,n,pos0_recv,pos1_recv,pos0,pos1
		real(f):: xmin_rem(dim),xmax_rem(dim),xi(dim),t1,t2,xmin_loc(dim),xmax_loc(dim),dr
		logical:: wait_for_phys_then_reloop
		integer:: ones1D(ntotal_loc+nphys_recv_all),halo_pindex_0(ntotal_loc+nphys_recv_all)
	
		! Initialization
		if ( allocated(halo_pindex) ) deallocate( halo_pindex )
		
		allocate( halo_pindex(ntotal_loc+nphys_recv_all,n_process_neighbour) )
			
		nhalo_send(:) = 0
	
		if (nphys_recv_all.eq.0) then
			wait_for_phys_then_reloop = .false.
		else
			wait_for_phys_then_reloop = .true.
		end if
	
		! Halo particle send location determination
		! First loop loops over currently held particles
		pos0 = 1
		pos1 = ntotal_loc
	
		! Begin search
		xmin_loc = bounds_glob(1:dim,procid+1) + scale_k*hsml
		xmax_loc = bounds_glob(dim+1:2*dim,procid+1) - scale_k*hsml
	1 	do i = pos0,pos1
			xi(:) = parts(i)%x(:)
			if ( any( [xi(:).le.xmin_loc(:) , xi(:).ge.xmax_loc(:)] ) ) then ! if particle is potentially neighbour's halo
				do j = 1,n_process_neighbour
					pid = proc_neighbour_list(j) + 1
					xmin_rem(:) = bounds_glob(1:dim,pid) - scale_k*hsml
					xmax_rem(:) = bounds_glob(dim+1:2*dim,pid) + scale_k*hsml
					if ( all( [xi(:).ge.xmin_rem(:) , xi(:).le.xmax_rem(:)] ) ) then
					
						nhalo_send(j) = nhalo_send(j) + 1
						halo_pindex(nhalo_send(j),j) = i
			
					end if
		
				end do
			end if
		end do
		
		! Waiting for physical particles to complete exchange if needed
		if (wait_for_phys_then_reloop) then
			pos0 = ntotal_loc+1
			pos1 = ntotal_loc+nphys_recv_all
			wait_for_phys_then_reloop = .false.
			call MPI_WAITALL(nrequest,request_in,status(:,1:nrequest),ierr) !wait for new physical particles to arrive
			goto 1
		end if
	
		! Posting non-blocking send for nhalo_send exchange
		do n = 1,n_process_neighbour
			pid = proc_neighbour_list(n)
			call MPI_IRECV(nhalo_recv(n),1,MPI_INTEGER,pid,0,MPI_COMM_WORLD,request_out(2*n-1),ierr)
			call MPI_ISEND(nhalo_send(n),1,MPI_INTEGER,pid,0,MPI_COMM_WORLD,request_out(2*n),ierr)
		end do
		
		! Creating indexed derived types to send halo particles to each neighbouring process
		ones1D(1:MAXVAL(nhalo_send)) = (/(1,i=1,MAXVAL(nhalo_send))/)
		do n = 1,n_process_neighbour
			if (nhalo_send(n) > 0) then
				halo_pindex_0(1:nhalo_send(n)) = halo_pindex(1:nhalo_send(n),n)-1
				! Type for initial exchange
				call MPI_TYPE_INDEXED(nhalo_send(n),ones1D,halo_pindex_0,halotype,halotype_indexed(n),ierr)
				call MPI_TYPE_COMMIT(halotype_indexed(n),ierr)
				! Type for subsequent exchanges
				call MPI_TYPE_INDEXED(nhalo_send(n),ones1D,halo_pindex_0,haloupdatetype,haloupdatetype_indexed(n),ierr)
				call MPI_TYPE_COMMIT(haloupdatetype_indexed(n),ierr)
			end if
		end do
		
		ntotal_loc = ntotal_loc + nphys_recv_all
		parts(1:ntotal_loc)%indloc = (/ (i, i=1,ntotal_loc) /) ! updating all the loca indices of new physical particles
	
		deallocate( PhysPackSend )
		
		! Wait for non-blocking send to complete
		call MPI_WAITALL(2*n_process_neighbour,request_out,status,ierr)
	
		! Non-blocking sends to exchange physical particles that have moved processes
		
		! stopping program if array bounds are exceeded
		nhalo_loc = SUM(nhalo_recv)
		if (ntotal_loc+nhalo_loc.gt.maxnloc) call error_msg(3,parts(i)%indloc)
		
		nrequest = 0
		! Creating communicator for indexed struct comms
		pos1_recv = ntotal_loc
		do n = 1,n_process_neighbour
	
			pid = proc_neighbour_list(n)
	
			if (nhalo_recv(n).gt.0) then

				pos0_recv = pos1_recv + 1
				pos1_recv = pos0_recv - 1 + nhalo_recv(n)
				
				nrequest = nrequest + 1
		
				call MPI_IRECV(parts(pos0_recv)%indglob,nhalo_recv(n),halotype,pid,0,MPI_COMM_WORLD,request_out(nrequest),ierr)
		
			end if
	
			if (nhalo_send(n).gt.0) then
			
				nrequest = nrequest + 1
		
				call MPI_ISEND(parts(1)%indglob,1,halotype_indexed(n),pid,0,MPI_COMM_WORLD,request_out(nrequest),ierr)
	
			end if
	
		end do
	
	end subroutine ORB_sendrecv_halo
	
	!==============================================================================================================================
	subroutine ORB_sendrecv_haloupdate(ki)
	! Reduced version of ORB_sendrecv_halo where no searching occurs (only the exchange)
	
		use globvar_para,	only: halotype_indexed,haloupdatetype_indexed,nhalo_recv,nhalo_send
		use param_para,		only: haloupdatetype
	
		implicit none
		integer,intent(in):: ki
		integer:: n,i,j,k,pos0_recv,pos1_recv,request(2*n_process_neighbour),pid,status(MPI_STATUS_SIZE,2*n_process_neighbour),n_request
	
		!3. halo particle send/receive ----------------------------------------------------------------------------------------------------
		n_request = 0
		pos1_recv = ntotal_loc
		do n = 1,n_process_neighbour
	
			pid = proc_neighbour_list(n)
	
			if (nhalo_recv(n).gt.0) then
				
				! Indices where receiving data will be placed
				pos0_recv = pos1_recv + 1
				pos1_recv = pos0_recv - 1 + nhalo_recv(n)
				
				n_request = n_request + 1
				call MPI_IRECV(parts(pos0_recv)%rho,nhalo_recv(n),haloupdatetype,pid,0,MPI_COMM_WORLD,request(n_request),ierr)
		
			end if
	
			if (nhalo_send(n).gt.0) then
			
				n_request = n_request + 1
				call MPI_ISEND(parts(1)%rho,1,haloupdatetype_indexed(n),pid,0,MPI_COMM_WORLD,request(n_request),ierr)
				
			end if
	
		end do
	
		! Freeing up halotype used just for first exchange
		if (ki.eq.2) then
			do i = 1,n_process_neighbour
				if (nhalo_send(i)>0) call MPI_TYPE_FREE(halotype_indexed(i),ierr)
			end do
		end if
	
		call MPI_WAITALL(n_request,request,status,ierr)
		
		! Freeing up halo type used for 2nd - 3rd exchanges
		if (ki.eq.4) then
			do i = 1,n_process_neighbour
				if (nhalo_send(i)>0) call MPI_TYPE_FREE(haloupdatetype_indexed(i),ierr)
			end do
		end if
		
	end subroutine ORB_sendrecv_haloupdate
	
end module ORB_sr_m
