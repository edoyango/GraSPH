!==========================================================
subroutine output( ) 
!==========================================================

	use globvar
	use param
	
	implicit none     
	integer:: i,d,n
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
      
end