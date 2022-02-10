module flink_list_m
	
	use datatypes, only: particles
	use globvar, only: parts,pairs,niac,ntotal,nvirt,maxinter,scale_k
	use param, only: dim,f,hsml
	
	use kernel_m, only: kernel
	
	type particleincellarray
		type(particles),pointer:: p
	end type particleincellarray
	
	public:: flink_list,particleincellarray
	private:: check_if_interact

contains
	subroutine flink_list( )
	! subroutine to search for particle interactions using cell-linked list.
		
		implicit none
		integer,parameter:: maxpcell=12
		integer:: ngridx(dim),i,j,k,d,icell,jcell,xi,yi
		real(f):: maxgridx(dim),mingridx(dim),dcell
		integer,allocatable:: pincell(:,:)
		type(particleincellarray),allocatable:: cells(:,:,:)
		integer,parameter:: sweep(2,4) = reshape((/ 1,-1,&
													1, 0,&
													0, 1,&
													1, 1/),(/2,4/))
	
		!Determining bounding box extents
		mingridx(:) = parts(1)%x(:)
		maxgridx(:) = parts(1)%x(:)
		do i = 2,ntotal+nvirt
			do d = 1,dim
				mingridx(d) = MIN(mingridx(d),parts(i)%x(d))
				maxgridx(d) = MAX(maxgridx(d),parts(i)%x(d))
			end do
		end do
		
		!Determining number of grid cells in each direction
		dcell = scale_k*hsml
		maxgridx(:) = maxgridx(:) + 2d0*dcell
		mingridx(:) = mingridx(:) - 2d0*dcell
		ngridx(:) = int((maxgridx(:) - mingridx(:))/dcell) + 1
		maxgridx(:) = mingridx(:) + ngridx(:)*dcell
		
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
	
	end subroutine flink_list
	
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
	
end module flink_list_m