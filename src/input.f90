!==========================================================
subroutine input(generate)
!==========================================================

	use globvar
	use param
	
	implicit none
	logical,intent(in):: generate
	integer:: i,j,d,n
	real(8):: xi,yi
	real(8),parameter:: dx=dxo,dy=dxo,xl=25d0,yl=25d0
	
	select case (generate)
	
		case (.false.)
			
			ntotal = mp*np
			
		case (.true.)
		
			n = 0
			do i = 1, mp
				do j = 1, np
					xi = i*dx - dx/2.d0
					yi = j*dy - dy/2.d0
					n = n + 1
					parts(n)%ind = n
					parts(n)%x(1) = xi
					parts(n)%x(2) = yi
					parts(n)%vx(:) = 0d0
					parts(n)%itype = 2
					parts(n)%rho = irho
					parts(n)%p = 0d0
				enddo
			enddo
			
			call write_ini_config
		
	end select

end subroutine input

!==================================================================================================
subroutine write_ini_config
! Subroutine for writing initial configuration data using intrinsic IO

	use globvar
	use param

	implicit none
	integer:: i,d,im

	open(10,file="C:\Users\edwar\Documents\outputdata/ini_xv.dat")
	open(20,file="C:\Users\edwar\Documents\outputdata/ini_state.dat")
	open(30,file="C:\Users\edwar\Documents\outputdata/ini_stress.dat")
	
	do i = 1, ntotal 
		write(10,1001) i, parts(i)%x(:), parts(i)%vx(:)
		write(20,1002) i, parts(i)%itype, hsml, mass, parts(i)%rho, parts(i)%p
	end do
	
	1001 format(1x, I5, 6(2x, e14.7)) 
	1002 format(1x, I5, 2x, I2, 8(2x, e14.7)) 
	1003 format(1x, I5, 2x, I2, 2x, e14.7) 
	
	close(1)
	close(2) 
	close(3) 
	
end subroutine write_ini_config