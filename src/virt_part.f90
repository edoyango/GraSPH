!=============================================================
subroutine virt_part(generate)
!=============================================================
	
	use globvar
	use param
	
	implicit none
	logical,intent(in):: generate
	integer:: i,k,d,n
	real(8),parameter:: dx = dxo, dy = dxo, xmin = 0d0, ymin = 0d0, xmax = 75d0, ymax = 40d0
	
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

end