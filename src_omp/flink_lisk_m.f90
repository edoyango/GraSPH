module flink_list_m
	
	use datatypes, only: particles
	use globvar, only: parts,pairs,niac,ntotal,nvirt,maxinter,scale_k
	use param, only: dim,f,hsml
	use omp_lib
	
	use kernel_m, only: kernel
	
	type particleincellarray
		type(particles),pointer:: p
	end type particleincellarray
	
	integer,private:: ngridx(dim),i,j,k,d,icell,jcell,kcell,xi,yi,zi,jth,tmpind
	real(f),private:: mingridx(dim),maxgridx(dim),dcell!,minx_loc(dim),maxx_loc(dim)
	
	public:: flink_list,particleincellarray
	private:: check_if_interact,flink_list2D,flink_list3D,bounding_box

contains
	
	!==============================================================================================================================
	subroutine flink_list
	! interface subroutine for 2D and 3D cases
		implicit none
		
		select case (dim)
			case(2)
				call flink_list2D
			case(3)
				call flink_list3D
		end select
		
	end subroutine flink_list
	
	!==============================================================================================================================
	subroutine flink_list2D( )
	! subroutine to search for particle interactions using cell-linked list.
		
		implicit none
		integer,parameter:: maxpcell=12
		integer,allocatable:: pincell(:,:)
		type(particleincellarray),allocatable:: cells(:,:,:)
		integer,parameter:: sweep(2,4) = reshape((/ 1,-1,&
													1, 0,&
													0, 1,&
													1, 1/),(/2,4/))
		real(f):: minx_loc(dim),maxx_loc(dim)
													
		call bounding_box(mingridx,maxgridx,ngridx,dcell)
		
		allocate( pincell(ngridx(1),ngridx(2)),&
					cells(maxpcell,ngridx(1),ngridx(2)) )
					
		!Mapping particles to grid cells
		pincell(:,:) = 0
		
		do i=1,ntotal+nvirt
			icell = int((parts(i)%x(1) - mingridx(1))/dcell) + 1
			jcell = int((parts(i)%x(2) - mingridx(2))/dcell) + 1 
			pincell(icell,jcell) = pincell(icell,jcell) + 1
			cells(pincell(icell,jcell),icell,jcell)%p => parts(i)
		enddo
		
		! Searching for interactions. Algorithm loops through every cell with indices icell,jcell.
		! First checks particles within cell are interacting, then checks particles interacting with
		! particles in adjacent cells. 
		niac = 0
		do icell = 2,ngridx(1)-2
			do jcell = 2,ngridx(2)-2
			
				! finding pairs within cell icell,jcell
				do i = 1,pincell(icell,jcell)-1
					do j = i+1,pincell(icell,jcell)
						call check_if_interact(cells(i,icell,jcell)%p,cells(j,icell,jcell)%p)
					end do
				end do
				
				! finding pairs between particles in cell icell,jcell and particles in cell xi,yi
				do k = 1,4
					xi = icell + sweep(1,k)
					yi = jcell + sweep(2,k)
					do i = 1,pincell(icell,jcell)
						do j = 1,pincell(xi,yi)
							call check_if_interact(cells(i,icell,jcell)%p,cells(j,xi,yi)%p)
						end do
					end do
				end do
				
			end do
		end do
	
	end subroutine flink_list2D
	
	!==============================================================================================================================
	subroutine flink_list3D( )
	! save as above, but for 3D
		
		implicit none
		integer,parameter:: maxpcell=64
		integer,allocatable:: pincell(:,:,:)
		type(particleincellarray),allocatable:: cells(:,:,:,:)
		real(f):: minx_loc(dim),maxx_loc(dim),t1
		integer,parameter:: sweep(3,13) = reshape((/-1,-1,-1,&
													-1,-1, 0,&
													-1,-1, 1,&
													-1, 0,-1,&
													-1, 0, 0,&
													-1, 0, 1,&
													-1, 1,-1,&
													-1, 1, 0,&
													-1, 1, 1,&
													 0,-1,-1,&
													 0,-1, 0,&
													 0,-1, 1,&
													 0, 0,-1/),(/3,13/))
		
		
		mingridx(:) = HUGE(1_f)
		maxgridx(:) = -HUGE(1_f)
		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(minx_loc,maxx_loc,t1)
		t1 = OMP_GET_WTIME()
		!Determining bounding box extents
		!call bounding_box(mingridx,maxgridx,ngridx,dcell)
		minx_loc(:) = HUGE(1_f)
		maxx_loc(:) = -HUGE(1_f)
		!$OMP DO COLLAPSE(2)
		do i = 1,ntotal+nvirt
			do d = 1,dim
				minx_loc(d) = MIN(minx_loc(d),parts(i)%x(d))
				maxx_loc(d) = MAX(maxx_loc(d),parts(i)%x(d))
			end do
		end do
		
		do d = 1,dim
			!$OMP ATOMIC
			mingridx(d) = MIN(mingridx(d),minx_loc(d))
			!$OMP ATOMIC
			maxgridx(d) = MAX(maxgridx(d),maxx_loc(d))
		end do
		!$OMP FLUSH(mingridx,maxgridx)
		
		!$OMP SINGLE
		dcell = scale_k*hsml
		mingridx(:) = mingridx(:) - 2_f*dcell
		maxgridx(:) = maxgridx(:) + 2_f*dcell
		ngridx(:) = int((maxgridx(:)-mingridx(:))/dcell) + 1
		maxgridx(:) = mingridx(:) + ngridx(:)*dcell
		!$OMP END SINGLE

		!$OMP FLUSH(ngridx)
		
		!$OMP SINGLE
		allocate( pincell(ngridx(1),ngridx(2),ngridx(3)),&
					cells(maxpcell,ngridx(1),ngridx(2),ngridx(3)) )
		!$OMP END SINGLE
		
		!$OMP DO COLLAPSE(3)
		do i = 1,ngridx(1)
			do j = 1,ngridx(2)
				do k = 1,ngridx(3)
					pincell(i,j,k) = 0
				end do
			end do
		end do
		!$OMP END DO
		
		!$OMP FLUSH(pincell,mingridx)
		
		!$OMP DO PRIVATE(icell,jcell,kcell,tmpind)
		do i=1,ntotal+nvirt
			icell = int((parts(i)%x(1) - mingridx(1))/dcell) + 1
			jcell = int((parts(i)%x(2) - mingridx(2))/dcell) + 1 
			kcell = int((parts(i)%x(3) - mingridx(3))/dcell) + 1 
			!$OMP ATOMIC CAPTURE
			pincell(icell,jcell,kcell) = pincell(icell,jcell,kcell) + 1
			tmpind = pincell(icell,jcell,kcell)
			!$OMP END ATOMIC
			cells(tmpind,icell,jcell,kcell)%p => parts(i)
		enddo
		!$OMP END DO
		write(*,*) OMP_GET_WTIME()-t1
		!$OMP SINGLE
		niac = 0
		do icell = 2,ngridx(1)-2
			do jcell = 2,ngridx(2)-2
				do kcell = 2,ngridx(3)-2
					if (pincell(icell,jcell,kcell) > 0) then
					
						! finding pairs within cell icell,jcell
						do i = 1,pincell(icell,jcell,kcell)-1
							do j = i+1,pincell(icell,jcell,kcell)
								call check_if_interact(cells(i,icell,jcell,kcell)%p,cells(j,icell,jcell,kcell)%p)
							end do
						end do
						
						! finding pairs between particles in cell icell,jcell and particles in cell xi,yi
						do k = 1,13
							xi = icell + sweep(1,k)
							yi = jcell + sweep(2,k)
							zi = kcell + sweep(3,k)
							do i = 1,pincell(icell,jcell,kcell)
								do j = 1,pincell(xi,yi,zi)
									call check_if_interact(cells(i,icell,jcell,kcell)%p,cells(j,xi,yi,zi)%p)
								end do
							end do
						end do
					
					endif
				end do
			end do
		end do
		!$OMP END SINGLE
		!$OMP END PARALLEL
	end subroutine flink_list3D
	
	!==============================================================================================================================
	subroutine check_if_interact(p_i,p_j)
	! subroutine to chekc if two particles are interacting and consequently adding to pair list
		
		implicit none
		type(particles),intent(in),target:: p_i,p_j
		real(f):: dxiac(dim),r
		
		! only consider interactions when real-real are involved
		if ( p_i%itype.eq.1 .or. p_j%itype.eq.1 ) then
			dxiac(:) = p_i%x(:) - p_j%x(:)
			r = SQRT(SUM(dxiac*dxiac))
			if (r < hsml*scale_k) then
				niac = niac + 1
				if (niac < maxinter) then
					pairs(niac)%i => p_i
					pairs(niac)%j => p_j
					call kernel(r,dxiac,hsml,pairs(niac)%w,pairs(niac)%dwdx(:))
				else
					print *,' >>> Error <<< : Too many interactions'
					stop
				end if
			end if
		end if
		
	end subroutine check_if_interact
	
	!==============================================================================================================================
	subroutine bounding_box(minx,maxx,ng,dc)
	
		implicit none
		real(f),intent(out):: minx(dim),maxx(dim),dc
		integer,intent(out):: ng(dim)
		real(f):: minx_loc(dim),maxx_loc(dim)
		minx_loc(:) = HUGE(1_f)
		maxx_loc(:) = -HUGE(1_f)
		!Determining bounding box extents
		!$OMP DO
		do i = 1,ntotal+nvirt
			do d = 1,dim
				minx_loc(d) = MIN(minx_loc(d),parts(i)%x(d))
				maxx_loc(d) = MAX(maxx_loc(d),parts(i)%x(d))
			end do
		end do
		!$OMP END DO
		!$OMP BARRIER
		write(*,*) minx_loc,maxx_loc
		!$OMP CRITICAL
		do d = 1,dim
			! !$OMP ATOMIC
			minx(d) = MIN(minx(d),minx_loc(d))
			! !$OMP ATOMIC
			maxx(d) = MAX(maxx(d),maxx_loc(d))
		end do
		!$OMP END CRITICAL
		
		!$OMP FLUSH(minx,maxx)
		
		!$OMP SINGLE
		!Determining number of grid cells in each direction
		dc = scale_k*hsml
		maxx(:) = maxx(:) + 2d0*dc
		minx(:) = minx(:) - 2d0*dc
		ng(:) = int((maxx(:) - minx(:))/dc) + 1
		maxx(:) = minx(:) + ng(:)*dc
		!$OMP END SINGLE
		
	end subroutine bounding_box
		
end module flink_list_m