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

end subroutine input

!==================================================================================================
subroutine write_ini_config
! Subroutine for writing initial configuration data using intrinsic IO

	use globvar
	use param

	implicit none
	integer:: i,d,im

	open(1,file=trim(output_directory)//"/ini_xv.dat")
	open(2,file=trim(output_directory)//"/ini_state.dat")
	open(3,file=trim(output_directory)//"/ini_stress.dat")
	
	do i = 1, ntotal 
		write(1,1001) i, parts(i)%itype, parts(i)%x(:), parts(i)%vx(:)
		write(2,1001) i, parts(i)%itype, parts(i)%rho, parts(i)%p
	end do
	
	1001 format(1x, I6, 1x, I2, 6(2x, e24.17))
	
	close(1)
	close(2) 
	close(3) 
	
end subroutine write_ini_config