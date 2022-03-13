module output_m
	
	use globvar, only: ntotal,nvirt,parts
	use param, only: output_directory
	
	integer,private:: i,d

contains

	!==============================================================================================================================
	subroutine output( ) 
	! Subroutine to write data to disk. Called intermittently, determined by save_step supplied at run-time
	
		use globvar,	only: itimestep,save_step
		use param,		only: output_phys,output_virt
		
		implicit none  
		integer:: n
		character(len=4)::number
		
		n = itimestep/save_step
		write(number,'(I4.4)') n
        
        open(1,file=trim(output_directory)//"/f_ind"//number//".dat",form='unformatted',access='stream',status='replace')
        do i = 1,ntotal
            write(1) parts(i)%ind,0,parts(i)%itype
        end do
		
		if (output_phys(1)) then
			
			open(1,file=trim(output_directory)//"/f_xv"//number//".dat",form='unformatted',access='stream',status='replace')
			
			do i = 1,ntotal
				write(1) parts(i)%x(:),parts(i)%vx(:)
			end do
			
			close(1)
			
		end if
		
		if (output_phys(2)) then
		
			open(2,file=trim(output_directory)//"/f_state"//number//".dat",form='unformatted',access='stream',status='replace')
			
			do i = 1,ntotal
				write(2) parts(i)%rho, parts(i)%p
			end do
			
			close(2)
			
		end if
        
        open(1,file=trim(output_directory)//"/v_ind"//number//".dat",form='unformatted',access='stream',status='replace')
        do i = ntotal+1,ntotal+nvirt
            write(1) parts(i)%ind,0,parts(i)%itype
        end do
		
		if (output_virt(1)) then
			
			open(3,file=trim(output_directory)//"/v_xv"//number//".dat",form='unformatted',access='stream',status='replace')
			
			do i = ntotal+1,ntotal+nvirt
				write(3) parts(i)%x(:), parts(i)%vx(:)
			end do
			
			close(3)
			
		end if
		
		if (output_virt(2)) then
			
			open(4,file=trim(output_directory)//"/v_state"//number//".dat",form='unformatted',access='stream',status='replace')
			
			do i = ntotal+1,ntotal+nvirt
				write(4) parts(i)%rho, parts(i)%p
			end do
			
			close(4)
			
		end if
		
	end subroutine output
	
	!==============================================================================================================================
	subroutine write_ini_config
	! Subroutine for writing initial configuration data using intrinsic IO
	
		implicit none
        
        open(1,file=trim(output_directory)//"/ini_ind.dat",form='unformatted',access='stream',status='replace')        
		open(2,file=trim(output_directory)//"/ini_xv.dat",form='unformatted',access='stream',status='replace')
		open(3,file=trim(output_directory)//"/ini_state.dat",form='unformatted',access='stream',status='replace')
		
		do i = 1, ntotal 
            write(1) parts(i)%ind, 0, parts(i)%itype
			write(2) parts(i)%x(:), parts(i)%vx(:)
			write(3) parts(i)%rho, parts(i)%p
		end do
		
		close(1); close(2); close(3)
		
	end subroutine write_ini_config
	
end module output_m
