module input_m
	
    use datatypes, only: particles
	use globvar, only: parts,pairs,ntotal,nvirt,nghos,scale_k,niac
	use param, only: dim,f,dxo,mp,np,op,pp,qp,rp,nlayer,irho,hsml,mass,rh0,gamma,c
	
	use output_m, only: write_ini_config
    
    real(f),parameter:: vxmin = 0._f, vymin = 0._f, vzmin = 0._f, &
        vxmax = vxmin + pp*dxo, vymax = vymin + qp*dxo, vzmax = vzmin + rp*dxo
    real(f),parameter:: rxmin = 0._f, rymin = 0._f, rzmin = 0._f, &
        rxmax = rxmin + mp*dxo, rymax = rymin + np*dxo, rzmax = rzmin + op*dxo
    
    integer,allocatable:: gind(:)
	
	public:: input,virt_part
	private:: input2D,input3D,virt_part2D,virt_part3D,vxmin,vymin,vzmin,vxmax,vymax,vzmax,rxmin,rymin,rzmin,rxmax,rymax,rzmax
	
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
		
		select case (generate)
		
			case (.false.)
				
				ntotal = mp*op
				
			case (.true.)
			
				n = 0
				do i = 1, mp
					do j = 1, op
						n = n + 1
						parts(n)%ind = n
						parts(n)%x(1) = (i-0.5_f)*dxo
						parts(n)%x(2) = (j-0.5_f)*dxo
						parts(n)%vx(:) = 0_f
						parts(n)%itype = 1
						parts(n)%rho = irho
						parts(n)%p = 0_f
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
        real(f),parameter:: K0=0.3_f/0.7_f
		
		select case (generate)
		
			case (.false.)
				
				ntotal = mp*np*op
				
			case (.true.)
			
				n = 0
                do k = 1,op
				
					do j = 1, np
						do i = 1, mp
							n = n + 1
							parts(n)%ind = n
							parts(n)%x(1) = (i-0.5_f)*dxo
							parts(n)%x(2) = (j-0.5_f)*dxo
							parts(n)%x(3) = (k-0.5_f)*dxo
							parts(n)%vx(:) = 0_f
							parts(n)%itype = 1
							parts(n)%rho = irho
							parts(n)%p = 0_f
!~                             parts(n)%sig(3) = -irho*9.81_f*(rzmax-parts(n)%x(3))
!~                             parts(n)%sig(1:2) = K0*parts(n)%sig(3)
                            parts(n)%sig(:) = 0._f
                            parts(n)%strain(:) = 0._f
                            parts(n)%pstrain(:) = 0._f
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
		
		select case (generate)
			
			case (.false.)
				
				nvirt = 2*rp + 2*pp
				
			case (.true.)
				
				n = ntotal
				
				!---Virtual particle on the lower boundary
				do i = 1, pp
					n = n + 1
					parts(n)%ind = n
					parts(n)%x(1) = vxmin + (i-0.5_f)*dxo
					parts(n)%x(2) = vymin - 0.5_f*dxo
					parts(n)%vx(:) = 0_f
				enddo
				
				!---Virtual particle on the upper boundary
				do i = 1, pp
					n = n + 1
					parts(n)%ind = n
					parts(n)%x(1) = vxmin + (i-0.5_f)*dxo
					parts(n)%x(2) = vymax - 1.5_f*dxo
					parts(n)%vx(:) = 0_f
				enddo
				
				!---Virtual particle on the left boundary
				do i = 1, rp
					n = n + 1
					parts(n)%ind = n
					parts(n)%x(1) = vxmin - 0.5_f*dxo
					parts(n)%x(2) = vymin + (i-1.5_f)*dxo
					parts(n)%vx(:) = 0_f
				enddo
				
				!---Virtual particle on the right boundary
				do i = 1, rp
					n = n + 1
					parts(n)%ind = n
					parts(n)%x(1) = vxmax + 0.5_f*dxo
					parts(n)%x(2) = vymin + (i-1.5_f)*dxo
					parts(n)%vx(:) = 0_f
				enddo
	
				parts(ntotal+1:ntotal+nvirt)%rho = irho
				parts(ntotal+1:ntotal+nvirt)%p = 0_f
				parts(ntotal+1:ntotal+nvirt)%itype = -1
				
		end select
	
	end subroutine virt_part2D
	
	!==============================================================================================================================
	subroutine virt_part3D(generate)
		
		implicit none
		logical,intent(in):: generate
		integer:: i,j,k,d,n
		
		select case (generate)
			
			case (.false.)
				
				nvirt = (pp+2*nlayer)*(qp+2*nlayer)*nlayer
				
			case (.true.)
				
				n = ntotal
				
				do i = 1-nlayer,pp+nlayer
                    do j = 1-nlayer,qp+nlayer
                        do k = 1,nlayer
                            n = n + 1
                            parts(n)%ind = n
                            parts(n)%x(1) = vxmin + (i-0.5_f)*dxo
                            parts(n)%x(2) = vymin + (j-0.5_f)*dxo
                            parts(n)%x(3) = vzmin - (k-0.5_f)*dxo
                            parts(n)%vx(:) = 0._f
                        end do
                    end do
                end do
                
				parts(ntotal+1:ntotal+nvirt)%rho = irho
				parts(ntotal+1:ntotal+nvirt)%p = 0_f
				parts(ntotal+1:ntotal+nvirt)%itype = -1
				
		end select
	
	end subroutine virt_part3D
    
    !==============================================================================================================================
    subroutine generate_ghost_part
        
        implicit none
        integer:: i,j,k,d,ig
        
        nghos = 0
        
        do i = 1,ntotal
            if (abs(parts(i)%x(1)-vxmin) < scale_k*hsml .and. parts(i)%x(1) > vxmin) then
                nghos = nghos + 1
                ig = ntotal+nvirt+nghos
                gind(nghos) = i
                parts(ig) = parts(i)
                parts(ig)%ind = ig
                parts(ig)%itype = 99
                parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmin
                parts(ig)%vx(1) = -parts(ig)%vx(1)
                parts(ig)%sig(4:6) = -parts(ig)%sig(4:6)
            end if
            if (abs(parts(i)%x(1)-vxmax) < scale_k*hsml .and. parts(i)%x(1) < vxmax) then
                nghos = nghos + 1
                ig = ntotal+nvirt+nghos
                gind(nghos) = i
                parts(ig) = parts(i)
                parts(ig)%ind = ig
                parts(ig)%itype = 99
                parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmax
                parts(ig)%vx(1) = -parts(ig)%vx(1)
                parts(ig)%sig(4:6) = -parts(ig)%sig(4:6)
            end if
            if (abs(parts(i)%x(2)-vymin) < scale_k*hsml .and. parts(i)%x(2) > vymin) then
                nghos = nghos + 1
                ig = ntotal+nvirt+nghos
                gind(nghos) = i
                parts(ig) = parts(i)
                parts(ig)%ind = ig
                parts(ig)%itype = 98
                parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymin
                parts(ig)%vx(2) = -parts(ig)%vx(2)
                parts(ig)%sig(4:6) = -parts(ig)%sig(4:6)
            end if
            if (abs(parts(i)%x(2)-vymax) < scale_k*hsml .and. parts(i)%x(2) < vymax) then
                nghos = nghos + 1
                ig = ntotal+nvirt+nghos
                gind(nghos) = i
                parts(ig) = parts(i)
                parts(ig)%ind = ig
                parts(ig)%itype = 98
                parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymax
                parts(ig)%vx(2) = -parts(ig)%vx(2)
                parts(ig)%sig(4:6) = -parts(ig)%sig(4:6)
            end if
            if ( (parts(i)%x(1)-vxmin)**2 + (parts(i)%x(2)-vymin)**2 < (scale_k*hsml)**2 .and. parts(i)%x(1) > vxmin) then
                nghos = nghos + 1
                ig = ntotal+nvirt+nghos
                gind(nghos) = i
                parts(ig) = parts(i)
                parts(ig)%ind = ig
                parts(ig)%itype = 97
                parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmin
                parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymin
                parts(ig)%vx(1) = -parts(ig)%vx(1)
                parts(ig)%vx(2) = -parts(ig)%vx(2)
                parts(ig)%sig(4:6) = -parts(ig)%sig(4:6)
            end if
            if ( (parts(i)%x(1)-vxmin)**2 + (parts(i)%x(2)-vymax)**2 < (scale_k*hsml)**2 .and. parts(i)%x(1) > vxmin) then
                nghos = nghos + 1
                ig = ntotal+nvirt+nghos
                gind(nghos) = i
                parts(ig) = parts(i)
                parts(ig)%ind = ig
                parts(ig)%itype = 97
                parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmin
                parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymax
                parts(ig)%vx(1) = -parts(ig)%vx(1)
                parts(ig)%vx(2) = -parts(ig)%vx(2)
                parts(ig)%sig(4:6) = -parts(ig)%sig(4:6)
            end if
            if ( (parts(i)%x(1)-vxmax)**2 + (parts(i)%x(2)-vymax)**2 < (scale_k*hsml)**2 .and. parts(i)%x(1) < vxmax) then
                nghos = nghos + 1
                ig = ntotal+nvirt+nghos
                gind(nghos) = i
                parts(ig) = parts(i)
                parts(ig)%ind = ig
                parts(ig)%itype = 97
                parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmax
                parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymax
                parts(ig)%vx(1) = -parts(ig)%vx(1)
                parts(ig)%vx(2) = -parts(ig)%vx(2)
                parts(ig)%sig(4:6) = -parts(ig)%sig(4:6)
            end if
            if ( (parts(i)%x(1)-vxmax)**2 + (parts(i)%x(2)-vymin)**2 < (scale_k*hsml)**2 .and. parts(i)%x(1) < vxmax) then
                nghos = nghos + 1
                ig = ntotal+nvirt+nghos
                gind(nghos) = i
                parts(ig) = parts(i)
                parts(ig)%ind = ig
                parts(ig)%itype = 97
                parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmax
                parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymin
                parts(ig)%vx(1) = -parts(ig)%vx(1)
                parts(ig)%vx(2) = -parts(ig)%vx(2)
                parts(ig)%sig(4:6) = -parts(ig)%sig(4:6)
            end if
        end do
        
    end subroutine generate_ghost_part
    
    !==============================================================================================================================
    subroutine update_ghost_part
        
        implicit none
        integer:: i,j,k,d,ig,ir
        
        do i = 1,nghos
            ig = ntotal+nvirt+i
            ir = gind(i)
            select case (parts(ig)%itype)
                case(99)
                    parts(ig)%rho = parts(ir)%rho
                    parts(ig)%p = parts(ir)%p
                    parts(ig)%vx(1) =-parts(ir)%vx(1)
                    parts(ig)%vx(2) = parts(ir)%vx(2)
                    parts(ig)%vx(3) = parts(ir)%vx(3)
                    parts(ig)%sig(1:dim) = parts(ir)%sig(1:dim)
                    parts(ig)%sig(4:6) = parts(ir)%sig(4:6)
                case(98)
                    parts(ig)%rho = parts(ir)%rho
                    parts(ig)%p = parts(ir)%p
                    parts(ig)%vx(1) = parts(ir)%vx(1)
                    parts(ig)%vx(2) =-parts(ir)%vx(2)
                    parts(ig)%vx(3) = parts(ir)%vx(3)
                    parts(ig)%sig(1:dim) = parts(ir)%sig(1:dim)
                    parts(ig)%sig(4:6) = parts(ir)%sig(4:6)
                case(97)
                    parts(ig)%rho = parts(ir)%rho
                    parts(ig)%p = parts(ir)%p
                    parts(ig)%vx(1) =-parts(ir)%vx(1)
                    parts(ig)%vx(2) =-parts(ir)%vx(2)
                    parts(ig)%vx(3) = parts(ir)%vx(3)
                    parts(ig)%sig(1:dim) = parts(ir)%sig(1:dim)
                    parts(ig)%sig(4:6) = parts(ir)%sig(4:6)
            end select
        end do
        
    end subroutine update_ghost_part
    
    !==============================================================================================================================
    subroutine update_virt_part
    
        implicit none
        integer:: i,j,k,d,iv,ir
        real(f):: tmp
        
        vw(:) = 0._f
        do i = 1,nvirt
            parts(ntotal+i)%rho = 0._f
            parts(ntotal+i)%vx(:) = 0._f
        end do
        
        do k = 1,niac
            
            if (pairs(k)%i%itype < 0 .and. pairs(k)%j%itype > 0) then
                tmp = mass*pairs(k)%w/pairs(k)%j%rho
                vw(pairs(k)%i%ind-ntotal) = vw(pairs(k)%i%ind-ntotal) + tmp
                pairs(k)%i%rho = pairs(k)%i%rho + mass*pairs(k)%w
                pairs(k)%i%vx(:) = pairs(k)%i%vx(:) - pairs(k)%j%vx(:)*tmp
            else if (pairs(k)%j%itype < 0 .and. pairs(k)%i%itype > 0) then
                tmp = mass*pairs(k)%w/pairs(k)%i%rho
                vw(pairs(k)%j%ind-ntotal) = vw(pairs(k)%j%ind-ntotal) + tmp
                pairs(k)%j%rho = pairs(k)%j%rho + mass*pairs(k)%w
                pairs(k)%j%vx(:) = pairs(k)%j%vx(:) - pairs(k)%i%vx(:)*tmp
            end if
        
        end do
        
        do i = 1,nvirt
            if (vw(i) .ne. 0._f) then
                parts(ntotal+i)%rho = parts(ntotal+i)%rho/vw(i)
                parts(ntotal+i)%vx(:) = parts(ntotal+i)%vx(:)/vw(i)
            else
                parts(ntotal+i)%rho = irho
                parts(ntotal+i)%vx(:) = 0._f
            end if
            parts(ntotal+i)%p = rh0*c**2*((parts(ntotal+i)%rho/rh0)**gamma-1_f)/gamma
        end do
        
    end subroutine update_virt_part
    
    !==============================================================================================================================
    subroutine virt_mirror(pr,pv)
    
        implicit none
        type(particles),intent(in):: pr
        type(particles),intent(inout):: pv
        real(f):: da,db,beta
        real(f),parameter:: beta_max = 5._f
        
        da = ABS(pr%x(3) - vzmin)
        db = ABS(pv%x(3) - vzmin)
        
        beta = MIN(1._f + db/da,beta_max)
        if (da==0._f) beta = beta_max
        
        pv%rho = pr%rho
        pv%p = pr%p
        pv%sig = pr%sig
        pv%vx(:) = (1._f-beta)*pr%vx(:)
        
    end subroutine virt_mirror

end module input_m
