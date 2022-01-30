!=================================================================================
subroutine con_density(ki,cdrhodt)
!=================================================================================

	use globvar
	use param
	
	implicit none
	integer,intent(in):: ki
	real(8),intent(out):: cdrhodt(ntotal+nvirt)
	integer:: i,j,k,d
	real(8):: dvx(dim),vcc
	
	cdrhodt(1:ntotal) = 0d0
	
	do k=1,niac      
	
		if (pairs(k)%i%itype.gt.0 .and. pairs(k)%j%itype.gt.0) then
		
			dvx(:) = pairs(k)%i%vx(:) - pairs(k)%j%vx(:) 
			
			vcc = DOT_PRODUCT(dvx(:),pairs(k)%dwdx(:))
			
			cdrhodt(pairs(k)%i%ind) = cdrhodt(pairs(k)%i%ind) + mass*vcc
			cdrhodt(pairs(k)%j%ind) = cdrhodt(pairs(k)%j%ind) + mass*vcc       
			if ((pairs(k)%i%ind==15505 .or. pairs(k)%j%ind==15505) .and. itimestep==1723) write(*,'(3(I6,1x),9(e14.7,1x))') ki,pairs(k)%i%ind,pairs(k)%j%ind,pairs(k)%i%vx,pairs(k)%j%vx,pairs(k)%dwdx
		endif
	
	enddo
	 
end subroutine con_density