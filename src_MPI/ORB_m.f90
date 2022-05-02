module ORB_m

	use globvar,		only: scale_k,parts,ntotal_loc
	use globvar_para,	only: procid,numprocs,ierr,MPI_ftype,repartition_mode,node_cax,node_cut
	use mpi
	use param,			only: f,dim,hsml
	
	public:: ORB
	private:: particle_grid,particle_grid2D,particle_grid3D,ORB_bounds,ORB_bounds2D,ORB_bounds3D,P_trim2D,P_trim3D,&
		data_to_next_node,ORB_boundingbox,potential_neighbour_process_search,subdomain_neighbour
	
contains	
	!==============================================================================================================================
	subroutine ORB
	! Container subroutine for the bulk of the ORB algorithm, including the initial exchange of physical and halo particles
		use globvar,		only: itimestep,ntotal,nhalo_loc,nvirt_loc,t_graph,t_dist
		use globvar_para,	only: n_reorients,n_parts,mintstep_bn_part,maxtstep_bn_part,mintstep_bn_reorient,maxtstep_bn_reorient,&
			box_ratio_previous,maxnode,bounds_glob,node_segment,prev_part_tstep,prev_reorient_tstep,nphys_send,nphys_recv,&
			nhalo_send,nhalo_recv,halotype_indexed,haloupdatetype_indexed,prev_load,n_process_neighbour,proc_neighbour_list,&
            VirtPack
		use param_para,		only: dcell_ORB,ORBcheck1,ORBcheck2,box_ratio_threshold
		use ORB_sr_m,		only: ORB_sendrecv_diffuse,ORB_sendrecv_halo
		use input_m,		only: virt_part
		
		implicit none
		real(f),parameter:: dcell=hsml*dcell_ORB
		integer:: d,i,j,k,ngridx(dim),nphys_recv_all,request_phys(2*numprocs),request_halo(2*numprocs),searchrange_ini(2),n_request,&
		status(MPI_STATUS_SIZE,4*numprocs),nphys_send_all,procrange_ini(2),tree_layers,gridind_ini(dim,2),diffusedepth
		real(f):: bounds_out(2*dim),mingridx_ini(dim),maxgridx_ini(dim),current_to_previous(dim,dim),box_ratio_current(dim,dim)
		
		!allocating partitioning arrays and initialising diagnostic variables -------------------------------------------------------------
		t_graph = t_graph - MPI_WTIME() ! commence timing of ORB algoirthm
		if (itimestep.eq.1) then
			tree_layers = CEILING(LOG(DBLE(numprocs))/LOG(2d0))
			maxnode = 2*2**tree_layers-1
			allocate( bounds_glob(2*dim,numprocs),&
					proc_neighbour_list(numprocs),&
					node_cax(maxnode),&
					node_segment(maxnode) )
		endif
		
		! Boundary Determiniation Algorithm ---------------------------------------------------------------------------------------
		repartition_mode = 1 !initially assumes no partition
		! only checks if boundary needs updating every 50-100 time-steps
		if ( (itimestep.eq.1) .or. ((itimestep-prev_part_tstep.ge.ORBcheck1) .and. (mod(itimestep-prev_part_tstep,ORBcheck2).eq.0)) ) then
		
			! checking if change in partilces on current process > 5%
			if (ntotal_loc.gt.prev_load+0.05_f*DBLE(ntotal)/DBLE(numprocs)) then
				repartition_mode = 2
			endif
		
			call MPI_ALLREDUCE(MPI_IN_PLACE,repartition_mode,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
		
			if ( (repartition_mode.gt.1) .or. (itimestep.eq.1) ) then
			
				n_parts = n_parts + 1
				if (itimestep.ne.1) then
					mintstep_bn_part = min(mintstep_bn_part,itimestep-prev_part_tstep)
					maxtstep_bn_part = max(maxtstep_bn_part,itimestep-prev_part_tstep)
				endif
				prev_part_tstep = itimestep
				
				! Creating particle-in-cell grid
				call particle_grid( ngridx,dcell,mingridx_ini,maxgridx_ini )
				
				! Calculating current aspect ratio.
				do d = 1,dim
					box_ratio_current(:,d) = DBLE(ngridx(d))/DBLE(ngridx(:))
				end do
				
				current_to_previous(:,:) = box_ratio_current(:,:)/box_ratio_previous(:,:)
				if ( any( current_to_previous(:,:).gt.1_f+box_ratio_threshold ) ) repartition_mode = 3
				
				!partition summary info
				if (repartition_mode.eq.3) then 
					box_ratio_previous(:,:) = box_ratio_current(:,:)
					if (itimestep.ne.1) then
						maxtstep_bn_reorient = max(maxtstep_bn_reorient,itimestep-prev_reorient_tstep)
						mintstep_bn_reorient = min(mintstep_bn_reorient,itimestep-prev_reorient_tstep)
					end if
					prev_reorient_tstep = itimestep
					n_reorients = n_reorients + 1
				endif
			
				! determine subdomain boundaries using particle distribution
				gridind_ini(:,1) = 1
				gridind_ini(:,2) = ngridx(:)
				procrange_ini(1) = 0
				procrange_ini(2) = numprocs-1
				bounds_out = ORB_bounds( gridind_ini,numprocs,1,procrange_ini,ntotal,dcell,mingridx_ini,maxgridx_ini )
	
				call subdomain_neighbour
                
                ! Local save subset of global virtual particles
                call virt_part(.true.)
				
				! Updating sizes of select arrays to account for potential changes in neighbour list size
				if (itimestep.ne.1) deallocate(nphys_send,nphys_recv,nhalo_send,nhalo_recv,halotype_indexed,haloupdatetype_indexed)
				allocate( nphys_send(n_process_neighbour),&
						nphys_recv(n_process_neighbour),&
						nhalo_send(n_process_neighbour),&
						nhalo_recv(n_process_neighbour),&
						halotype_indexed(n_process_neighbour),&
						haloupdatetype_indexed(n_process_neighbour) )
			
			endif
		
		endif
		t_graph = t_graph + MPI_WTIME( ) ! conclude timing of ORB algorithm
	
		! Particle distribution (physical, halo) ----------------------------------------------------------------------------------
		t_dist = t_dist - MPI_WTIME() ! commence timing of particle distribution
		
		! physical particle distribution
		diffusedepth = 0
		searchrange_ini(:) = (/1,ntotal_loc/)
		i = ORB_sendrecv_diffuse( diffusedepth,searchrange_ini,n_request,request_phys,nphys_recv_all )
		
		! halo particle distribution
		call ORB_sendrecv_halo( request_phys,request_halo,nphys_recv_all,n_request )
		
		! update virtual particles
!~ 		call virt_part(.true.)
        parts(ntotal_loc+nhalo_loc+1:ntotal_loc+nhalo_loc+nvirt_loc) = VirtPack(1:nvirt_loc)
		
		if (repartition_mode.gt.1) prev_load = ntotal_loc
		
		do i = ntotal_loc+1,ntotal_loc+nhalo_loc
			parts(i)%indloc = i
		end do
		
		! wait for halo particle distribution to complete
		call MPI_WAITALL(n_request,request_halo(1:n_request),status(:,1:n_request),ierr)
		
		parts(ntotal_loc+1:ntotal_loc+nhalo_loc)%itype = 2
		t_dist = t_dist + MPI_WTIME()
		
	end subroutine ORB
	
	!==============================================================================================================================
	subroutine particle_grid(ngridx,dcell,mingridx,maxgridx)
	! interface subroutine
		
		implicit none
		real(f),intent(in):: dcell
		integer,intent(out):: ngridx(dim)
		real(f),intent(out):: mingridx(dim),maxgridx(dim)
		
		select case (dim)
			case (2)
				call particle_grid2D(ngridx,dcell,mingridx,maxgridx)
			case (3)
				call particle_grid3D(ngridx,dcell,mingridx,maxgridx)
		end select
		
	end subroutine particle_grid
	
	!==============================================================================================================================
	function ORB_bounds( gridind_in,nprocs_in,node_in,procrange_in,ntotal_in,dcell,mingridx_in,maxgridx_in ) result(bounds_out)
	! interface subroutine
		
		implicit none
		integer,intent(in):: gridind_in(dim,2),nprocs_in,node_in,procrange_in(2),ntotal_in
		real(f),intent(in):: dcell,mingridx_in(dim),maxgridx_in(dim)
		real(f):: bounds_out(2*dim)
		
		select case (dim)
			case (2)
				bounds_out = ORB_bounds2D( gridind_in,nprocs_in,node_in,procrange_in,ntotal_in,dcell,mingridx_in,maxgridx_in )
			case (3)
				bounds_out = ORB_bounds3D( gridind_in,nprocs_in,node_in,procrange_in,ntotal_in,dcell,mingridx_in,maxgridx_in )
		end select
		
	end function ORB_bounds
	
	!==============================================================================================================================
	subroutine particle_grid2D( ngridx,dcell,mingridx,maxgridx )
	! Subroutine to create a uniform rectangular grid with square cells, and counting the number of particles contained within each
	! cell. Each MPI process holds a global copy of the entire grid, in preperation for ORB
		
		use globvar_para,	only: pincell_ORB_2D
		
		implicit none
		real(f),intent(in):: dcell
		integer,intent(out):: ngridx(2)
		real(f),intent(out):: mingridx(2),maxgridx(2)
		integer:: i,j,k,d,icell,jcell,n_nonzerocells,n_nonzerocells_perprocess(numprocs),n_nonzerocells_total,pid,displ(numprocs),&
			cellmins(2),cellmaxs(2),cellrange(2),request,status(MPI_STATUS_SIZE)
		integer:: sendcount,recvcount(numprocs),Plist_size
		real(f):: maxx(2),minx(2)
		integer,allocatable:: Plist_loc(:,:),Plist_all(:,:)
		
		call ORB_boundingbox(ngridx,minx,maxx,mingridx,maxgridx,cellmins,cellmaxs,cellrange,dcell,2)
		
		! Creating list of non-zero grid cells ------------------------------------------------------------------------------------
		Plist_size = MIN(ntotal_loc,PRODUCT(cellrange(:)))
		allocate( pincell_ORB_2D(ngridx(1),ngridx(2)),&
				Plist_loc(3,Plist_size) )
				
		pincell_ORB_2D(cellmins(1):cellmaxs(1),cellmins(2):cellmaxs(2)) = 0
		
		! Counting particles in each cell
		n_nonzerocells = 0
		do i = 1,ntotal_loc
			icell = int((parts(i)%x(1)-mingridx(1))/dcell) + 1
			jcell = int((parts(i)%x(2)-mingridx(2))/dcell) + 1
			if ( pincell_ORB_2D(icell,jcell).eq.0) then
				n_nonzerocells = n_nonzerocells + 1
				Plist_loc(1,n_nonzerocells) = icell
				Plist_loc(2,n_nonzerocells) = jcell
			endif
			pincell_ORB_2D(icell,jcell) = pincell_ORB_2D(icell,jcell) + 1
		enddo
		
		! Exchanging info on how many non-zero cells each process has
		call MPI_IALLGATHER(n_nonzerocells,1,MPI_INTEGER,n_nonzerocells_perprocess,1,MPI_INTEGER,MPI_COMM_WORLD,request,ierr)
		
		! Populating Plist_loc with non-zero cell entries
		do i = 1,n_nonzerocells
			Plist_loc(3,i) = Pincell_ORB_2D(Plist_loc(1,i),Plist_loc(2,i))
		enddo
		
		call MPI_WAIT(request,status,ierr)
		
		n_nonzerocells_total = 0
		do pid = 1,numprocs
			displ(pid) = n_nonzerocells_total
			n_nonzerocells_total = n_nonzerocells_total + n_nonzerocells_perprocess(pid)
		enddo
		
		allocate( Plist_all(3,n_nonzerocells_total) )
		
		! Collecting all process' Plist_loc
		displ = 3*displ
		sendcount = 3*n_nonzerocells
		recvcount = 3*n_nonzerocells_perprocess(:)
		call MPI_IALLGATHERV(Plist_loc,sendcount,MPI_INTEGER,Plist_all,recvcount,displ,MPI_INTEGER,MPI_COMM_WORLD,request,ierr)
		
		pincell_ORB_2D(:,:) = 0
		
		call MPI_WAIT(request,status,ierr)
		
		! populating the number of particles per cell
		do i = 1,n_nonzerocells_total
			pincell_ORB_2D(Plist_all(1,i),Plist_all(2,i)) = pincell_ORB_2D(Plist_all(1,i),Plist_all(2,i)) + Plist_all(3,i)
		enddo
		
		deallocate( Plist_loc,Plist_all )
	
	end subroutine particle_grid2D
	
	!==============================================================================================================================
	recursive function ORB_bounds2D( gridind_in,nprocs_in,node_in,procrange_in,ntotal_in,dcell,mingridx_in,maxgridx_in ) &
		result(bounds_out)
	! Recursive function that performs the 'bisection' part of the ORB algorithm
	
		use globvar_para,	only: leaf_node,bounds_glob,pincell_ORB_2D
		use param_para,		only: bound_extend
	
		implicit none
		integer,intent(in):: gridind_in(dim,2),node_in,nprocs_in,procrange_in(2),ntotal_in
		real(f),intent(in):: mingridx_in(dim),maxgridx_in(dim),dcell
		integer:: i,j,k,d,node_out,gridind_out(dim,2),nprocs_out,ntotal_out,procrange_out(2),n_p,cax,np_per_node,pincol,&
			ngridx_trim(dim)
		real(f):: t1,t2,bounds_out(2*dim)

		!determining cut axis. 1 = x, 2 = y ---------------------------------------------------------------------------------------
		if (repartition_mode.eq.3) then
			call P_trim2D(gridind_in,ngridx_trim)
			cax = 1
			if (ngridx_trim(2) .gt. ngridx_trim(1)) cax = 2
			node_cax(node_in) = cax
		else
			cax = node_cax(node_in)
		endif

		!cut location -------------------------------------------------------------------------------------------------------------
		i = gridind_in(cax,1) - 1
		n_p = 0
		np_per_node = ceiling(real(nprocs_in)/2)/real(nprocs_in)*ntotal_in
		do while (n_p .lt. np_per_node)
			i = i + 1
			pincol = 0
			if (cax .eq. 1) then
				do j = gridind_in(2,1),gridind_in(2,2)
					pincol = pincol + pincell_ORB_2D(i,j)
				enddo
			elseif (cax.eq.2) then
				do j = gridind_in(1,1),gridind_in(1,2)
					pincol = pincol + pincell_ORB_2D(j,i)
				enddo
			endif
			n_p = n_p + pincol
		enddo
	
		if ( (np_per_node-(n_p-pincol).lt.n_p-np_per_node) ) then
			i = i - 1
			n_p = n_p - pincol
		endif

		!saving output information ------------------------------------------------------------------------------------------------
		call data_to_next_node(procrange_in,gridind_in,nprocs_in,procid,n_p,i,node_in,ntotal_in,cax,procrange_out,gridind_out,&
		nprocs_out,ntotal_out,node_out)
		
		!travelling to child node/saving boundary information ---------------------------------------------------------------------
		if (nprocs_out .ne. 1) then
			bounds_out = ORB_bounds2D( gridind_out,nprocs_out,node_out,procrange_out,ntotal_out,dcell,mingridx_in,maxgridx_in )
		else
			leaf_node = node_out
			
			bounds_out(1:dim) = mingridx_in(:) + (gridind_out(:,1)-1)*dcell
			bounds_out(dim+1:2*dim) = mingridx_in(:) + gridind_out(:,2)*dcell
			
			if (repartition_mode.eq.3) call potential_neighbour_process_search(leaf_node)
			
			if (node_cut(1).eq.0) bounds_out(3) = bounds_out(3) + bound_extend*scale_k*hsml
			if (node_cut(2).eq.0) bounds_out(1) = bounds_out(1) - bound_extend*scale_k*hsml
			if (node_cut(3).eq.0) bounds_out(4) = bounds_out(4) + bound_extend*scale_k*hsml
			if (node_cut(4).eq.0) bounds_out(2) = bounds_out(2) - bound_extend*scale_k*hsml
			
			call MPI_ALLGATHER(bounds_out,2*dim,MPI_ftype,bounds_glob,2*dim,MPI_ftype,MPI_COMM_WORLD,ierr)
			
			deallocate( pincell_ORB_2D )
		
		endif
  
	end function ORB_bounds2D
	
	!==============================================================================================================================
	subroutine P_trim2D(gridind_in,ngridx_trim)
	! Trims particle-in-cell grid so as to obtain minimal bounding boxes to obtain accurate cut axis orientations
		
		use globvar_para,	only: pincell_ORB_2D
		
		implicit none
		integer,intent(in):: gridind_in(dim,2)
		integer,intent(out):: ngridx_trim(dim)
		integer:: i,j,k,d,newi(2),newj(2),oldi(2),oldj(2)
		
		oldi(:) = gridind_in(1,:)
		oldj(:) = gridind_in(2,:)
	
		!Trimming input grid ------------------------------------------------------------------------------------------------------
		!reducing x-length of grid
		!finding new start x-index
		newstarti: do i = oldi(1),oldi(2)
			if ( any(pincell_ORB_2D(i,oldj(1):oldj(2)).ne.0) ) then
				newi(1) = i
				exit newstarti
			end if
		end do newstarti

		!finding new end x-index
 		newendi: do i = oldi(2),newi(1),-1
			if ( any(pincell_ORB_2D(i,oldj(1):oldj(2)).ne.0) ) then
				newi(2) = i
				exit newendi
			end if
		end do newendi

		!reducing y-length of grid
		!finding new start y-index
 		newstartj: do j = oldj(1),oldj(2)
			if ( any(pincell_ORB_2D(newi(1):newi(2),j).ne.0) ) then
				newj(1) = j
				exit newstartj
			end if
		end do newstartj

		!finding new end y-index
 		newendj: do j = oldj(2),newj(1),-1
			if ( any(pincell_ORB_2D(newi(1):newi(2),j).ne.0) ) then
				newj(2) = j
				exit newendj
			end if
		end do newendj

 		ngridx_trim(1) = newi(2) - newi(1) + 1
		ngridx_trim(2) = newj(2) - newj(1) + 1

	end subroutine P_trim2D
	
	!==============================================================================================================================
	subroutine particle_grid3D( ngridx,dcell,mingridx,maxgridx )
	! Subroutine to create a uniform rectangular grid with square cells, and counting the number of particles contained within each
	! cell. Each MPI process holds a global copy of the entire grid, in preperation for ORB
		
		use globvar_para,	only: pincell_ORB_3D
		
		implicit none
		real(f),intent(in):: dcell
		integer,intent(out):: ngridx(:)
		real(f),intent(out):: mingridx(:),maxgridx(:)
		integer:: i,j,k,d,n,icell,jcell,kcell,n_nonzerocells,n_nonzerocells_perprocess(numprocs),n_nonzerocells_total,pid,&
				displ(numprocs),cellmins(3),cellmaxs(3),cellrange(3),request,status(MPI_STATUS_SIZE)
		real(f):: minx(3),maxx(3)
		integer:: sendcount,recvcount(numprocs),Plist_size
		integer,allocatable:: Plist_loc(:,:),Plist_all(:,:)
		
		call ORB_boundingbox(ngridx,minx,maxx,mingridx,maxgridx,cellmins,cellmaxs,cellrange,dcell,3)
		
		! Creating list of non-zero grid cells ------------------------------------------------------------------------------------
		Plist_size = MIN(ntotal_loc,PRODUCT(cellrange(:)))
		allocate( pincell_ORB_3D(ngridx(1),ngridx(2),ngridx(3)),&
				Plist_loc(4,Plist_size) )
				
		pincell_ORB_3D(cellmins(1):cellmaxs(1),cellmins(2):cellmaxs(2),cellmins(3):cellmaxs(3)) = 0
		
		! Counting particles in each cell
		n_nonzerocells = 0
		do i = 1,ntotal_loc
			icell = int((parts(i)%x(1)-mingridx(1))/dcell) + 1
			jcell = int((parts(i)%x(2)-mingridx(2))/dcell) + 1
			kcell = int((parts(i)%x(3)-mingridx(3))/dcell) + 1
			if ( pincell_ORB_3D(icell,jcell,kcell).eq.0) then
				n_nonzerocells = n_nonzerocells + 1
				Plist_loc(1,n_nonzerocells) = icell
				Plist_loc(2,n_nonzerocells) = jcell
				Plist_loc(3,n_nonzerocells) = kcell
			endif
			pincell_ORB_3D(icell,jcell,kcell) = pincell_ORB_3D(icell,jcell,kcell) + 1
		enddo
		
		! Exchanging info on how many non-zero cells each process has
		call MPI_IALLGATHER(n_nonzerocells,1,MPI_INTEGER,n_nonzerocells_perprocess,1,MPI_INTEGER,MPI_COMM_WORLD,request,ierr)
		
		! Populating Plist_loc with non-zero cell entries
		do i = 1,n_nonzerocells
			Plist_loc(4,i) = Pincell_ORB_3D(Plist_loc(1,i),Plist_loc(2,i),Plist_loc(3,i))
		enddo
		
		call MPI_WAIT(request,status,ierr)
		
		n_nonzerocells_total = 0
		do pid = 1,numprocs
			displ(pid) = n_nonzerocells_total
			n_nonzerocells_total = n_nonzerocells_total + n_nonzerocells_perprocess(pid)
		enddo
		
		allocate( Plist_all(4,n_nonzerocells_total) )
		
		! Collecting all process' Plist_loc
		displ = 4*displ
		sendcount = 4*n_nonzerocells
		recvcount = 4*n_nonzerocells_perprocess(:)
		call MPI_IALLGATHERV(Plist_loc,sendcount,MPI_INTEGER,Plist_all,recvcount,displ,MPI_INTEGER,MPI_COMM_WORLD,request,ierr)
		
		pincell_ORB_3D(:,:,:) = 0
		
		call MPI_WAIT(request,status,ierr)
		
		! populating the number of particles per cell
		do i = 1,n_nonzerocells_total
			pincell_ORB_3D(Plist_all(1,i),Plist_all(2,i),Plist_all(3,i)) = &
				pincell_ORB_3D(Plist_all(1,i),Plist_all(2,i),Plist_all(3,i)) + Plist_all(4,i)
		enddo
		
		deallocate( Plist_loc,Plist_all )
	
	end subroutine particle_grid3D
	
	!==============================================================================================================================
	recursive function ORB_bounds3D( gridind_in,nprocs_in,node_in,procrange_in,ntotal_in,dcell,mingridx_in,maxgridx_in ) &
		result(bounds_out)
	! Recursive function that performs the 'bisection' part of the ORB algorithm
	
		use globvar_para,	only: leaf_node,bounds_glob,pincell_ORB_3D
		use param_para,		only: bound_extend
	
		implicit none
		integer,intent(in):: gridind_in(dim,2),node_in,nprocs_in,procrange_in(2),ntotal_in
		real(f),intent(in):: mingridx_in(dim),maxgridx_in(dim),dcell
		integer:: i,j,k,d,node_out,gridind_out(dim,2),nprocs_out,ntotal_out,procrange_out(2),n_p,cax,np_per_node,pincol,&
		ngridx_trim(dim),A(3)
		real(f):: t1,t2,bounds_out(2*dim)

		!determining cut axis. 1 = x, 2 = y ---------------------------------------------------------------------------------------
		if (repartition_mode.eq.3) then
			call P_trim3D(gridind_in,ngridx_trim)
			A(1) = ngridx_trim(2)*ngridx_trim(3); A(2) = ngridx_trim(1)*ngridx_trim(3); A(3) = ngridx_trim(1)*ngridx_trim(2)
			if ( (A(1).le.A(2)) .and. (A(1).le.A(3)) ) cax = 1
			if ( (A(2).le.A(1)) .and. (A(2).le.A(3)) ) cax = 2
			if ( (A(3).le.A(1)) .and. (A(3).le.A(2)) ) cax = 3
			node_cax(node_in) = cax
		else
			cax = node_cax(node_in)
		endif

		!cut location -------------------------------------------------------------------------------------------------------------
		i = gridind_in(cax,1) - 1
		n_p = 0
		np_per_node = ceiling(real(nprocs_in)/2)/real(nprocs_in)*ntotal_in
		do while (n_p .lt. np_per_node)
			i = i + 1
			pincol = 0
			if (cax .eq. 1) then
				pincol = SUM(pincell_ORB_3D(&
					i,&
					gridind_in(2,1):gridind_in(2,2),&
					gridind_in(3,1):gridind_in(3,2)))
			elseif (cax .eq. 2) then
				pincol = SUM(pincell_ORB_3D(&
					gridind_in(1,1):gridind_in(1,2),&
					i,&
					gridind_in(3,1):gridind_in(3,2)))
			else
				pincol = SUM(pincell_ORB_3D(&
					gridind_in(1,1):gridind_in(1,2),&
					gridind_in(2,1):gridind_in(2,2),&
					i))
			endif
			n_p = n_p + pincol
		enddo
	
		if ( (np_per_node-(n_p-pincol).lt.n_p-np_per_node) ) then
			i = i - 1
			n_p = n_p - pincol
		endif
		
		!saving output information ------------------------------------------------------------------------------------------------
		call data_to_next_node(procrange_in,gridind_in,nprocs_in,procid,n_p,i,node_in,ntotal_in,cax,procrange_out,gridind_out,&
		nprocs_out,ntotal_out,node_out)
		
		!travelling to child node/saving boundary information ---------------------------------------------------------------------
		if (nprocs_out .ne. 1) then
			bounds_out = ORB_bounds3D( gridind_out,nprocs_out,node_out,procrange_out,ntotal_out,dcell,mingridx_in,maxgridx_in )
		else
			leaf_node = node_out
			
			bounds_out(1:dim) = mingridx_in(:) + (gridind_out(:,1)-1)*dcell
			bounds_out(dim+1:2*dim) = mingridx_in(:) + gridind_out(:,2)*dcell
			
			if (repartition_mode.eq.3) call potential_neighbour_process_search(leaf_node)
			
			if (node_cut(1).eq.0) bounds_out(4) = bounds_out(4) + bound_extend*scale_k*hsml
			if (node_cut(2).eq.0) bounds_out(1) = bounds_out(1) - bound_extend*scale_k*hsml
			if (node_cut(3).eq.0) bounds_out(5) = bounds_out(5) + bound_extend*scale_k*hsml
			if (node_cut(4).eq.0) bounds_out(2) = bounds_out(2) - bound_extend*scale_k*hsml
			if (node_cut(5).eq.0) bounds_out(6) = bounds_out(6) + bound_extend*scale_k*hsml
			if (node_cut(6).eq.0) bounds_out(3) = bounds_out(3) - bound_extend*scale_k*hsml

			call MPI_ALLGATHER(bounds_out,2*dim,MPI_ftype,bounds_glob,2*dim,MPI_ftype,MPI_COMM_WORLD,ierr)
			
			deallocate( pincell_ORB_3D )
		
		endif
  
	end function ORB_bounds3D
	
	!==============================================================================================================================
	subroutine P_trim3D(gridind_in,ngridx_trim)
	! Trims particle-in-cell grid so as to obtain minimal bounding boxes to obtain accurate cut axis orientations
		
		use globvar_para,	only: pincell_ORB_3D
		
		implicit none
		integer,intent(in):: gridind_in(dim,2)
		integer,intent(out):: ngridx_trim(dim)
		integer:: i,j,k,d,newi(2),newj(2),newk(2),oldi(2),oldj(2),oldk(2)
		
		oldi(:) = gridind_in(1,:)
		oldj(:) = gridind_in(2,:)
		oldk(:) = gridind_in(3,:)
	
		!Trimming input grid ------------------------------------------------------------------------------------------------------
		!reducing x-length of grid
		!finding new start x-index
		newstarti: do i = oldi(1),oldi(2)
			if ( any(pincell_ORB_3D(i,oldj(1):oldj(2),oldk(1):oldk(2)).ne.0) ) then
				newi(1) = i
				exit newstarti
			end if
		end do newstarti

		!finding new end x-index
 		newendi: do i = oldi(2),newi(1),-1
			if ( any(pincell_ORB_3D(i,oldj(1):oldj(2),oldk(1):oldk(2)).ne.0) ) then
				newi(2) = i
				exit newendi
			end if
		end do newendi

		!reducing y-length of grid
		!finding new start y-index
 		newstartj: do j = oldj(1),oldj(2)
			if ( any(pincell_ORB_3D(newi(1):newi(2),j,oldk(1):oldk(2)).ne.0) ) then
				newj(1) = j
				exit newstartj
			end if
		end do newstartj

		!finding new end y-index
 		newendj: do j = oldj(2),newj(1),-1
			if ( any(pincell_ORB_3D(newi(1):newi(2),j,oldk(1):oldk(2)).ne.0) ) then
				newj(2) = j
				exit newendj
			end if
		end do newendj
		
		!reducing z-length of grid
		!finding new start z-index
 		newstartk: do k = oldk(1),oldk(2)
			if ( any(pincell_ORB_3D(newi(1):newi(2),oldj(1):oldj(2),k).ne.0) ) then
				newk(1) = k
				exit newstartk
			end if
		end do newstartk

		!finding new end z-index
 		newendk: do k = oldk(2),newk(1),-1
			if ( any(pincell_ORB_3D(newi(1):newi(2),oldj(1):oldj(2),k).ne.0) ) then
				newk(2) = k
				exit newendk
			end if
		end do newendk

 		ngridx_trim(1) = newi(2) - newi(1) + 1
		ngridx_trim(2) = newj(2) - newj(1) + 1
		ngridx_trim(3) = newk(2) - newk(1) + 1

	end subroutine P_trim3D

	!==============================================================================================================================
	subroutine data_to_next_node(procrange_in,gridind_in,nprocs_in,procid,n_p,ind,node_in,ntotal_in,cax,procrange_out,gridind_out,&
	nprocs_out,ntotal_out,node_out)
	! helper subroutine for ORB_bounds to determine which node the process is travelling to
	
		implicit none
		integer,intent(in):: procrange_in(2),gridind_in(dim,2),nprocs_in,procid,n_p,ind,node_in,ntotal_in,cax
		integer,intent(out):: procrange_out(2),gridind_out(dim,2),nprocs_out,ntotal_out,node_out
		integer:: procrange_lo(2),procrange_hi(2)
		
		procrange_lo(1) = procrange_in(1)
		procrange_lo(2) = procrange_in(1) + ceiling(real(nprocs_in)/2)-1
		procrange_hi(1) = procrange_lo(2) + 1
		procrange_hi(2) = procrange_lo(1) + nprocs_in - 1
		gridind_out(:,1) = gridind_in(:,1)
		gridind_out(:,2) = gridind_in(:,2)
		
		if (procid.le.procrange_lo(2)) then
			node_out = 2*node_in
			procrange_out = procrange_lo
			ntotal_out = n_p
			gridind_out(cax,1) = gridind_in(cax,1)
			gridind_out(cax,2) = ind
		else
			node_out = 2*node_in + 1
			procrange_out = procrange_hi
			ntotal_out = ntotal_in - n_p
			gridind_out(cax,1) = ind + 1
			gridind_out(cax,2) = gridind_in(cax,2)
		endif
		
		nprocs_out = procrange_out(2) - procrange_out(1) + 1
		
	end subroutine data_to_next_node
	
	!==============================================================================================================================
	subroutine ORB_boundingbox(ngridx,minx,maxx,mingridx,maxgridx,cellmins,cellmaxs,cellrange,dcell,ndims)
	! shared subroutine for 2D, 3D particle_grid subroutines to determine bounding box of all particles.
		
		integer,intent(in):: ndims
		real(f),intent(in):: dcell
		integer,intent(out):: ngridx(ndims),cellmins(ndims),cellmaxs(ndims),cellrange(ndims)
		real(f),intent(out):: minx(ndims),maxx(ndims),mingridx(ndims),maxgridx(ndims)
		integer:: i,d,request(2),status(MPI_STATUS_SIZE,2)
		
		! Local max, min, in each direction ---------------------------------------------------------------------------------------
		minx(1:ndims) = parts(1)%x(1:ndims)
		maxx(1:ndims) = parts(1)%x(1:ndims)
		do i = 2,ntotal_loc
			do d = 1,ndims
				if (parts(i)%x(d).gt.maxx(d)) maxx(d) = parts(i)%x(d)
				if (parts(i)%x(d).lt.minx(d)) minx(d) = parts(i)%x(d)
			end do
		end do
		
		! Global max, min, in each direction --------------------------------------------------------------------------------------
		call MPI_IALLREDUCE(maxx(:),maxgridx(:),ndims,MPI_ftype,MPI_MAX,MPI_COMM_WORLD,request(1),ierr) !finding max over all processes
		call MPI_IALLREDUCE(minx(:),mingridx(:),ndims,MPI_ftype,MPI_MIN,MPI_COMM_WORLD,request(2),ierr) !finding min over all processes
		
		call MPI_WAITALL(2,request,status,ierr)
		
		mingridx(:) = mingridx(:) - dcell
		maxgridx(:) = maxgridx(:) + dcell
		
		! Number of grid cells and adjusting max extent in each direction ---------------------------------------------------------
		ngridx(:) = int((maxgridx(:)-mingridx(:))/dcell) + 1
		maxgridx(:) = mingridx(:) + dcell*ngridx(:)
		
		! Reducing search area by locating indices in which local particles are contained within
		cellmins(:) = int((minx(:) - mingridx(:))/dcell) + 1
		cellmaxs(:) = int((maxx(:) - mingridx(:))/dcell) + 1
		cellrange(:) = cellmaxs(:) - cellmins(:) + 1
		
	end subroutine ORB_boundingbox
	
	!==============================================================================================================================
	subroutine potential_neighbour_process_search(ID_node)
	! Used to determine which nodes are located at the edge of the domain in each direction.
	! IGNORE BELOW
	! Finds any potential neighbours by exploiting the binary tree data structure
	! First finds nodes in the tree that define the edge of the local process. Then finds processes that also use that node as an 
	! edge but on the opposite side. E.g. if local process has node 2 as an east edge, potential east neighbours of the local 
	! process are the processes that use node 2 as a west edge.
		
		implicit none
		integer:: i,j,k,d,pid,ID_node,start_node,blah,direction_skip
		
		!Initialization
		if (.not.allocated(node_cut)) allocate( node_cut(2*dim) )
		
		!finding nodes that define edge for current process		  
		node_cut(:) = 0
		do while (ID_node.ne.1)
		
			if (mod(ID_node,2).eq.0) then
				ID_node = ID_node/2
				if ( (node_cax(ID_node).eq.1) .and. (node_cut(1).eq.0) ) node_cut(1) = ID_node !right face of process defined by 'parent'
				if ( (node_cax(ID_node).eq.2) .and. (node_cut(3).eq.0) ) node_cut(3) = ID_node !top face of process defined by 'parent'
				if ( (node_cax(ID_node).eq.3) .and. (node_cut(5).eq.0) ) node_cut(5) = ID_node !top face of process defined by 'parent'
			else
				ID_node = (ID_node-1)/2
				if ( (node_cax(ID_node).eq.1) .and. (node_cut(2).eq.0) ) node_cut(2) = ID_node !left of process defined by 'parent'
				if ( (node_cax(ID_node).eq.2) .and. (node_cut(4).eq.0) ) node_cut(4) = ID_node !bottom of process defined by 'parent'
				if ( (node_cax(ID_node).eq.3) .and. (node_cut(6).eq.0) ) node_cut(6) = ID_node !bottom face of process defined by 'parent'
			endif
		
		enddo
	
	end subroutine potential_neighbour_process_search
	
	!==============================================================================================================================
	subroutine subdomain_neighbour
	!creates list of adjacent subdomains for the local subdomain by searching potential neighbours and seeing if they overlap
		
		use globvar_para,	only: proc_neighbour_list,n_process_neighbour,bounds_glob
		
		implicit none
		integer:: i,j,k,d,pid
		real(f):: bounds_loc_min(dim),bounds_loc_max(dim),bounds_rem_min(dim),bounds_rem_max(dim)
		
		!Checking if potential neighbours overlap with local process. Overlap = adjacent.
		n_process_neighbour = 0
		bounds_loc_min(1:dim) = bounds_glob(1:dim,procid+1) - 0.5_f*scale_k*hsml
		bounds_loc_max(1:dim) = bounds_glob(dim+1:2*dim,procid+1) + 0.5_f*scale_k*hsml
		neighboursearch: do pid = 1,numprocs
			if (pid.ne.procid+1) then
				bounds_rem_min(:) = bounds_glob(1:dim,pid) - 0.5_f*scale_k*hsml
				bounds_rem_max(:) = bounds_glob(dim+1:2*dim,pid) + 0.5_f*scale_k*hsml
				! if local and remote process' extended boundaries don't overlap, check next process
				if ( any([bounds_rem_max(:) < bounds_loc_min(:),bounds_rem_min(:) > bounds_loc_max(:)]) ) cycle neighboursearch
				n_process_neighbour = n_process_neighbour + 1
				proc_neighbour_list(n_process_neighbour) = pid - 1
			endif
	 	end do neighboursearch
	
	end subroutine subdomain_neighbour
	
end module ORB_m
