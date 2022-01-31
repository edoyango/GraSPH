module output_m
	
	use globvar
	use param
	
	integer:: i,d

contains

	!==============================================================================================================================
	subroutine output( ) 
	! Subroutine to write data to disk. Called intermittently, determined by save_step supplied at run-time
	
		use globvar
		use param
		
		implicit none  
		integer:: n
		character(len=4)::number
		
		n = itimestep/save_step
		write(number,1000) n
		
		if (output_phys(1)) then
			
			open(1,file=trim(output_directory)//"/f_xv"//number//".dat")
			
			do i = 1,ntotal
				write(1,1001) i, parts(i)%itype, parts(i)%x(:), parts(i)%vx(:)
			end do
			
			close(1)
			
		end if
		
		if (output_phys(2)) then
		
			open(2,file=trim(output_directory)//"/f_state"//number//".dat")
			
			do i = 1,ntotal
				write(2,1001) i, parts(i)%itype, parts(i)%rho, parts(i)%p
			end do
			
			close(2)
			
		end if
		
		if (output_virt(1)) then
			
			open(3,file=trim(output_directory)//"v_xv"//number//".dat")
			
			do i = ntotal+1,ntotal+nvirt
				write(3,1001) i, parts(i)%itype, parts(i)%x(:), parts(i)%vx(:)
			end do
			
			close(3)
			
		end if
		
		if (output_virt(2)) then
			
			open(4,file=trim(output_directory)//"v_state"//number//".dat")
			
			do i = ntotal+1,ntotal+nvirt
				write(4,1001) i, parts(i)%itype, parts(i)%rho, parts(i)%p
			end do
			
			close(4)
			
		end if
		
		1000 format(I4.4)
		1001 format(1x, I6, 1x, I2, 6(2x, e14.7))
		
	end subroutine output
	
	!==============================================================================================================================
	subroutine write_ini_config
	! Subroutine for writing initial configuration data using intrinsic IO
	
		implicit none
	
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
	
end module output_m