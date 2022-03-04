module input_m
	
	use globvar, only: parts,ntotal,nvirt
	use param, only: dim,f,dxo,mp,np,op,pp,qp,rp,irho
	
	use output_m, only: write_ini_config
	
	public:: input,virt_part
	private:: input2D,input3D,virt_part2D,virt_part3D
	
contains
	
	!==============================================================================================================================
	subroutine input(generate)
	! interface subroutine for initial particle config in the 2D and 3D case
	
		implicit none
		logical,intent(in):: generate
	
		select case (dim)
			case(2)
				call input2D(generate)
			case(3)
				call input3D(generate)
		end select
	
	end subroutine input
	
	!==============================================================================================================================
	subroutine virt_part(generate)
	! interface subroutine for virtual particle config in the 2D and 3D case
	
		implicit none
		logical,intent(in):: generate
		
		select case (dim)
			case(2)
				call virt_part2D(generate)
			case(3)
				call virt_part3D(generate)
		end select
		
	end subroutine virt_part

	!==============================================================================================================================
	subroutine input2D(generate)
	
		implicit none
		logical,intent(in):: generate
		integer:: i,j,d,n
		real(f):: xi,yi
		real(f),parameter:: dx=dxo,dy=dxo,xl=25d0,yl=25d0
		
		select case (generate)
		
			case (.false.)
				
				ntotal = mp*np
				
			case (.true.)
			
				n = 0
				do i = 1, mp
					do j = 1, np
						n = n + 1
						parts(n)%ind = n
						parts(n)%x(1) = (i-0.5d0)*dx
						parts(n)%x(2) = (j-0.5d0)*dy
						parts(n)%vx(:) = 0d0
						parts(n)%itype = 1
						parts(n)%rho = irho
						parts(n)%p = 0d0
					enddo
				enddo
				
				call write_ini_config
			
		end select
	
	end subroutine input2D
	
	!==============================================================================================================================
	subroutine input3D(generate)
	
		implicit none
		logical,intent(in):: generate
		integer:: i,j,k,d,n
		real(f):: xi,yi,zi
		real(f),parameter:: dx=dxo,dy=dxo,dz=dxo,xl=25d0,yl=2d0,zl=40d0
		
		select case (generate)
		
			case (.false.)
				
				ntotal = mp*np*op
				
			case (.true.)
			
				n = 0
				do i = 1, mp
					do j = 1, np
						do k = 1,op
							n = n + 1
							parts(n)%ind = n
							parts(n)%x(1) = (i-0.5d0)*dx
							parts(n)%x(2) = (j-0.5d0)*dy
							parts(n)%x(3) = (k-0.5d0)*dz
							parts(n)%vx(:) = 0d0
							parts(n)%itype = 1
							parts(n)%rho = irho
							parts(n)%p = 0d0
						end do
					end do
				end do
				
				call write_ini_config
			
		end select
	
	end subroutine input3D
	
	!==============================================================================================================================
	subroutine virt_part2D(generate)
		
		implicit none
		logical,intent(in):: generate
		integer:: i,j,k,d,n
		real(f),parameter:: dx = dxo, dy = dxo, xmin = 0d0, ymin = 0d0, xmax = 75d0, ymax = 40d0
		
		select case (generate)
			
			case (.false.)
				
				nvirt = 2*op + 2*pp
				
			case (.true.)
				
				n = ntotal
				
				!---Virtual particle on the lower boundary
				do i = 1, op
					n = n + 1
					parts(n)%ind = n
					parts(n)%x(1) = xmin + (i-0.5d0)*dx
					parts(n)%x(2) = ymin - 0.5d0*dy
					parts(n)%vx(:) = 0d0
				enddo
				
				!---Virtual particle on the upper boundary
				do i = 1, op
					n = n + 1
					parts(n)%ind = n
					parts(n)%x(1) = xmin + (i-0.5d0)*dx
					parts(n)%x(2) = ymax - 1.5d0*dy
					parts(n)%vx(:) = 0d0
				enddo
				
				!---Virtual particle on the left boundary
				do i = 1, pp
					n = n + 1
					parts(n)%ind = n
					parts(n)%x(1) = xmin - 0.5d0*dx
					parts(n)%x(2) = ymin + (i-1.5d0)*dy
					parts(n)%vx(:) = 0d0
				enddo
				
				!---Virtual particle on the right boundary
				do i = 1, pp
					n = n + 1
					parts(n)%ind = n
					parts(n)%x(1) = xmax + 0.5d0*dx
					parts(n)%x(2) = ymin + (i-1.5d0)*dy
					parts(n)%vx(:) = 0d0
				enddo
	
				parts(ntotal+1:ntotal+nvirt)%rho = irho
				parts(ntotal+1:ntotal+nvirt)%p = 0d0
				parts(ntotal+1:ntotal+nvirt)%itype = -1
				
		end select
	
	end subroutine virt_part2D
	
	!==============================================================================================================================
	subroutine virt_part3D(generate)
		
		implicit none
		logical,intent(in):: generate
		integer:: i,j,k,d,n
		real(f),parameter:: dx = dxo, dy = dxo, dz = dxo, xmin = 0d0, ymin = 0d0, zmin = 0d0,&
		xmax = xmin + pp*dx, ymax = ymin + qp*dy, zmax = zmin + rp*dz
		
		select case (generate)
			
			case (.false.)
				
				nvirt = pp*qp*2 + rp*qp*2 + pp*rp*2 + 4*rp
				
			case (.true.)
				
				n = ntotal
				
				!---Virtual particle on the bottom face
				do i = 1, pp
					do j = 1,qp
						n = n + 1
						parts(n)%ind = n
						parts(n)%x(1) = xmin + (i-0.5d0)*dx
						parts(n)%x(2) = ymin + (j-0.5d0)*dy
						parts(n)%x(3) = zmin - 0.5d0*dz
						parts(n)%vx(:) = 0d0
					end do
				end do
				
				!---Virtual particle on top face
				do i = 1,pp
					do j = 1,qp
						n = n + 1
						parts(n)%ind = n
						parts(n)%x(1) = xmin + (i-0.5d0)*dx
						parts(n)%x(2) = ymin + (j-0.5d0)*dy
						parts(n)%x(3) = zmax - 1.5d0*dz
						parts(n)%vx(:) = 0d0
					end do
				end do
				
				!---Virtual particle on the west face
				do j = 1,qp
					do k = 1,rp
						n = n + 1
						parts(n)%ind = n
						parts(n)%x(1) = xmin - 0.5d0*dx
						parts(n)%x(2) = ymin + (j-0.5d0)*dy
						parts(n)%x(3) = zmin + (k-1.5d0)*dz
						parts(n)%vx(:) = 0d0
					end do
				end do
				
				!---Virtual particle on the east face
				do j = 1,qp
					do k = 1,rp
						n = n + 1
						parts(n)%ind = n
						parts(n)%x(1) = xmax + 0.5d0*dx
						parts(n)%x(2) = ymin + (j-0.5d0)*dy
						parts(n)%x(3) = zmin + (k-1.5d0)*dz
						parts(n)%vx(:) = 0d0
					end do
				end do
				
				!---Virtual particle on the south face
				do i = 0,pp+1
					do k = 1,rp
						n = n + 1
						parts(n)%ind = n
						parts(n)%x(1) = xmin + (i-0.5d0)*dx
						parts(n)%x(2) = ymin - 0.5d0*dy
						parts(n)%x(3) = zmin + (k-1.5d0)*dz
						parts(n)%vx(:) = 0d0
					end do
				end do
				
				!---Virtual particle on the north face
				do i = 0,pp+1
					do k = 1,rp
						n = n + 1
						parts(n)%ind = n
						parts(n)%x(1) = xmin + (i-0.5d0)*dx
						parts(n)%x(2) = ymax + 0.5d0*dy
						parts(n)%x(3) = zmin + (k-1.5d0)*dz
						parts(n)%vx(:) = 0d0
					end do
				end do
	
				parts(ntotal+1:ntotal+nvirt)%rho = irho
				parts(ntotal+1:ntotal+nvirt)%p = 0d0
				parts(ntotal+1:ntotal+nvirt)%itype = -1
				
		end select
	
	end subroutine virt_part3D

end module input_m