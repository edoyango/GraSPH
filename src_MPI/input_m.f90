module input_m
	
	use globvar,		only: ntotal,nvirt,ntotal_loc,nhalo_loc,nvirt_loc,parts,scale_k,maxnloc
	use globvar_para,	only: procid,numprocs,bounds_glob
	use mpi
	use param,			only: dim,irho,dxo,f,hsml,mp,np,op,pp,qp,rp
	use error_msg_m,	only: error_msg
	use output_m,		only: write_ini_config
	
	public:: input,virt_part
	
contains
	
	!==============================================================================================================================
	subroutine input(generate)
	! Generates initial physical particle configuration.
	! 2 cases: return only number of particles retrieved, or generating the particles
	
		implicit none
		logical,intent(in):: generate
		integer:: i,j,k,d,n,n_loc,n_loc_i,n_start,n_done
		real(f):: xi,yi
		real(f),parameter:: xl = 25_f, yl = 25_f
	
		select case (generate)
			
			case (.false.)
			
				ntotal = mp*np*op
			
			case (.true.)
				
				! how many particles to generate per process
				n_loc_i = ceiling(dble(ntotal)/numprocs)
				if (procid.eq.numprocs-1) then
					n_loc = ntotal - (numprocs-1)*n_loc_i
				else
					n_loc = n_loc_i
				end if
				n_start = procid*n_loc_i + 1
				n_done = n_start + n_loc_i - 1
				
				! stopping program if array bounds are exceeded
				if ( (procid.eq.0) .and. (n_loc.gt.maxnloc) ) call error_msg(1,1)
		
				! intitial setup
				n = 0
				ntotal_loc = 0
				do i = 1,mp
					do j = 1,np
						do k = 1,op
							n = n + 1 ! tracking total number of particles generated
							! Only generating particles assigned to process
							if ( (n.ge.n_start) .and. (n.le.n_done) ) then
								ntotal_loc = ntotal_loc + 1
								parts(ntotal_loc)%indglob = n
								parts(ntotal_loc)%indloc = ntotal_loc
								parts(ntotal_loc)%x(1) = (i-0.5_f)*dxo
								parts(ntotal_loc)%x(2) = (j-0.5_f)*dxo
								parts(ntotal_loc)%x(3) = (k-0.5_f)*dxo
								parts(ntotal_loc)%vx(:) = 0_f
								parts(ntotal_loc)%itype = 1
								parts(ntotal_loc)%rho = irho
								parts(ntotal_loc)%p = 0_f
							end if
						end do
					end do
				end do
				
				call write_ini_config
			
		end select
		
	end subroutine input
	
	!==============================================================================================================================
	subroutine virt_part(generate)
	! Generates the virtual particle configuration. Can change over time or remain static
	! 2 cases: return only number of particles retrieved, or generating the particles
		
		implicit none
		integer:: i,j,k,d,n
		real(f):: xi(dim),xmin_loc(dim),xmax_loc(dim)
		logical,intent(in):: generate
		real(f),parameter:: dx = dxo, dy = dxo, dz = dxo, xmin = 0_f, ymin = 0_f, zmin = 0_f, xmax = xmin + dxo*pp, &
		ymax = ymin + dxo*qp, zmax = zmin + dxo*rp
		
		select case (generate)
			
			case (.false.)
				
				nvirt = pp*qp*2 + rp*qp*2 + pp*rp*2 + 4*rp
				
			case (.true.)
				
				xmin_loc(:) = bounds_glob(1:dim,procid+1) - scale_k*hsml
				xmax_loc(:) = bounds_glob(dim+1:2*dim,procid+1) + scale_k*hsml
				
				nvirt_loc = 0
				n = ntotal ! counter used to track particle indices
				
				!---Virtual particle on the bottom face
				do i = 1, pp
					do j = 1,qp
						n = n + 1
						xi(1) = xmin + (i-0.5_f)*dx
						xi(2) = ymin + (j-0.5_f)*dy
						xi(3) = zmin - 0.5_f*dz
						if ( all( [xi(:).ge.xmin_loc(:),xi(:).le.xmax_loc(:)] ) ) then
							nvirt_loc = nvirt_loc + 1
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%indglob = n
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%indloc = ntotal_loc+nhalo_loc+nvirt_loc
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%itype = -1
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%x(:) = xi(:)
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%vx(:) = 0_f
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%rho = irho
						end if
					end do
				end do
				
				!---Virtual particle on top face
				do i = 1,pp
					do j = 1,qp
						n = n + 1
						xi(1) = xmin + (i-0.5_f)*dx
						xi(2) = ymin + (j-0.5_f)*dy
						xi(3) = zmax - 1.5_f*dz
						if ( all( [xi(:).ge.xmin_loc(:),xi(:).le.xmax_loc(:)] ) ) then
							nvirt_loc = nvirt_loc + 1
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%indglob = n
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%indloc = ntotal_loc+nhalo_loc+nvirt_loc
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%itype = -1
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%x(:) = xi(:)
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%vx(:) = 0_f
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%rho = irho
						end if
					end do
				end do
				
				!---Virtual particle on the west face
				do j = 1,qp
					do k = 1,rp
						n = n + 1
						xi(1) = xmin - 0.5_f*dx
						xi(2) = ymin + (j-0.5_f)*dy
						xi(3) = zmin + (k-1.5_f)*dz
						if ( all( [xi(:).ge.xmin_loc(:),xi(:).le.xmax_loc(:)] ) ) then
							nvirt_loc = nvirt_loc + 1
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%indglob = n
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%indloc = ntotal_loc+nhalo_loc+nvirt_loc
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%itype = -1
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%x(:) = xi(:)
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%vx(:) = 0_f
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%rho = irho
						end if
					end do
				end do
				
				!---Virtual particle on the east face
				do j = 1,qp
					do k = 1,rp
						n = n + 1
						xi(1) = xmax + 0.5_f*dx
						xi(2) = ymin + (j-0.5_f)*dy
						xi(3) = zmin + (k-1.5_f)*dz
						if ( all( [xi(:).ge.xmin_loc(:),xi(:).le.xmax_loc(:)] ) ) then
							nvirt_loc = nvirt_loc + 1
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%indglob = n
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%indloc = ntotal_loc+nhalo_loc+nvirt_loc
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%itype = -1
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%x(:) = xi(:)
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%vx(:) = 0_f
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%rho = irho
						end if
					end do
				end do
				
				!---Virtual particle on the south face
				do i = 0,pp+1
					do k = 1,rp
						n = n + 1
						xi(1) = xmin + (i-0.5_f)*dx
						xi(2) = ymin - 0.5_f*dy
						xi(3) = zmin + (k-1.5_f)*dz
						if ( all( [xi(:).ge.xmin_loc(:),xi(:).le.xmax_loc(:)] ) ) then
							nvirt_loc = nvirt_loc + 1
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%indglob = n
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%indloc = ntotal_loc+nhalo_loc+nvirt_loc
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%itype = -1
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%x(:) = xi(:)
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%vx(:) = 0_f
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%rho = irho
						end if
					end do
				end do
				
				!---Virtual particle on the north face
				do i = 0,pp+1
					do k = 1,rp
						n = n + 1
						xi(1) = xmin + (i-0.5_f)*dx
						xi(2) = ymax + 0.5_f*dy
						xi(3) = zmin + (k-1.5_f)*dz
						if ( all( [xi(:).ge.xmin_loc(:),xi(:).le.xmax_loc(:)] ) ) then
							nvirt_loc = nvirt_loc + 1
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%indglob = n
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%indloc = ntotal_loc+nhalo_loc+nvirt_loc
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%itype = -1
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%x(:) = xi(:)
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%vx(:) = 0_f
							parts(ntotal_loc+nhalo_loc+nvirt_loc)%rho = irho
						end if
					end do
				end do
				
		end select
						
	end subroutine virt_part

end module input_m
