module input_m
	
    use datatypes,      only: particles
	use globvar,		only: ntotal,nvirt,ntotal_loc,nhalo_loc,nvirt_loc,nghos_loc,parts,scale_k,maxnloc
	use globvar_para,	only: procid,numprocs,bounds_glob
	use mpi
	use param,			only: dim,irho,dxo,f,hsml,mp,np,op,pp,qp,rp,nlayer,rh0,c,gamma,f,g
	use error_msg_m,	only: error_msg
	use output_m,		only: write_ini_config
    
    real(f),parameter:: vxmin = 0._f, vymin = 0._f, vzmin = 0._f, &
        vxmax = vxmin + pp*dxo, vymax = vymin + qp*dxo, vzmax = vzmin + rp*dxo
    real(f),parameter:: rxmin = 0._f, rymin = 0._f, rzmin = 0._f, &
        rxmax = rxmin + mp*dxo, rymax = rymin + np*dxo, rzmax = rzmin + op*dxo
        
    integer,allocatable:: gind(:)
	
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
								parts(ntotal_loc)%x(1) = rxmin + (i-0.5_f)*dxo
								parts(ntotal_loc)%x(2) = rymin + (j-0.5_f)*dxo
								parts(ntotal_loc)%x(3) = rzmin + (k-0.5_f)*dxo
								parts(ntotal_loc)%vx(:) = 0._f
								parts(ntotal_loc)%itype = 1
                                parts(ntotal_loc)%p = irho*g*(rzmax-parts(ntotal_loc)%x(3))
								parts(ntotal_loc)%rho = rh0*(parts(ntotal_loc)%p*gamma/(rh0*c**2)+1._f)**(1._f/gamma)
							end if
						end do
					end do
				end do
				
		end select
		
	end subroutine input
	
	!==============================================================================================================================
	subroutine virt_part(generate)
	! Generates the virtual particle configuration. Can change over time or remain static
	! 2 cases: return only number of particles retrieved, or generating the particles
		
		implicit none
		integer:: i,j,k,d,n,n_loc,n_loc_i,n_start,n_done
		logical,intent(in):: generate
		
		select case (generate)
			
			case (.false.)
				
				nvirt = nlayer*(2*nlayer+pp)*(2*nlayer+qp)
				
			case (.true.)
            
                ! how many particles to generate per process
				n_loc_i = ceiling(dble(nvirt)/numprocs)
				if (procid.eq.numprocs-1) then
					n_loc = nvirt - (numprocs-1)*n_loc_i
				else
					n_loc = n_loc_i
				end if
				n_start = procid*n_loc_i + 1
				n_done = n_start + n_loc_i - 1
				
				nvirt_loc = 0
				n = ntotal ! counter used to track particle indices
				
				!---Virtual particle on the bottom face
				do i = 1-nlayer, pp+nlayer
					do j = 1-nlayer,qp+nlayer
                        do k = 1,nlayer
                            n = n + 1
                            if ( (n-ntotal.ge.n_start) .and. (n-ntotal.le.n_done) ) then
                                nvirt_loc = nvirt_loc + 1
                                parts(ntotal_loc+nvirt_loc)%indglob = n
                                parts(ntotal_loc+nvirt_loc)%indloc = ntotal_loc+nvirt_loc
                                parts(ntotal_loc+nvirt_loc)%itype = -1
                                parts(ntotal_loc+nvirt_loc)%x(1) = vxmin + (i-0.5_f)*dxo
                                parts(ntotal_loc+nvirt_loc)%x(2) = vymin + (j-0.5_f)*dxo
                                parts(ntotal_loc+nvirt_loc)%x(3) = vzmin - (k-0.5_f)*dxo
                                parts(ntotal_loc+nvirt_loc)%vx(:) = 0._f
                                parts(ntotal_loc+nvirt_loc)%rho = irho
                            end if
                        end do
					end do
				end do

		end select

	end subroutine virt_part
    
    !==============================================================================================================================
    subroutine generate_ghost_part
        
        implicit none
        integer:: i,j,k,d,ig
        
        nghos_loc = 0
        
        do i = 1,ntotal_loc+nvirt_loc+nhalo_loc
            if (parts(i)%itype > 0) then
                if (abs(parts(i)%x(1)-vxmin) < scale_k*hsml .and. parts(i)%x(1) > vxmin) then
                    nghos_loc = nghos_loc + 1
                    ig = ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc
                    gind(nghos_loc) = i
                    parts(ig) = parts(i)
                    parts(ig)%indloc = ig
                    parts(ig)%itype = 99
                    parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmin
                    parts(ig)%vx(1) = -parts(ig)%vx(1)
                end if
                if (abs(parts(i)%x(1)-vxmax) < scale_k*hsml .and. parts(i)%x(1) < vxmax) then
                    nghos_loc = nghos_loc + 1
                    ig = ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc
                    gind(nghos_loc) = i
                    parts(ig) = parts(i)
                    parts(ig)%indloc = ig
                    parts(ig)%itype = 99
                    parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmax
                    parts(ig)%vx(1) = -parts(ig)%vx(1)
                end if
                if (abs(parts(i)%x(2)-vymin) < scale_k*hsml .and. parts(i)%x(2) > vymin) then
                    nghos_loc = nghos_loc + 1
                    ig = ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc
                    gind(nghos_loc) = i
                    parts(ig) = parts(i)
                    parts(ig)%indloc = ig
                    parts(ig)%itype = 98
                    parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymin
                    parts(ig)%vx(2) = -parts(ig)%vx(2)
                end if
                if (abs(parts(i)%x(2)-vymax) < scale_k*hsml .and. parts(i)%x(2) < vymax) then
                    nghos_loc = nghos_loc + 1
                    ig = ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc
                    gind(nghos_loc) = i
                    parts(ig) = parts(i)
                    parts(ig)%indloc = ig
                    parts(ig)%itype = 98
                    parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymax
                    parts(ig)%vx(2) = -parts(ig)%vx(2)
                end if
                if ( (parts(i)%x(1)-vxmin)**2 + (parts(i)%x(2)-vymin)**2 < (scale_k*hsml)**2 .and. parts(i)%x(1) > vxmin) then
                    nghos_loc = nghos_loc + 1
                    ig = ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc
                    gind(nghos_loc) = i
                    parts(ig) = parts(i)
                    parts(ig)%indloc = ig
                    parts(ig)%itype = 97
                    parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmin
                    parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymin
                    parts(ig)%vx(1) = -parts(ig)%vx(1)
                    parts(ig)%vx(2) = -parts(ig)%vx(2)
                end if
                if ( (parts(i)%x(1)-vxmin)**2 + (parts(i)%x(2)-vymax)**2 < (scale_k*hsml)**2 .and. parts(i)%x(1) > vxmin) then
                    nghos_loc = nghos_loc + 1
                    ig = ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc
                    gind(nghos_loc) = i
                    parts(ig) = parts(i)
                    parts(ig)%indloc = ig
                    parts(ig)%itype = 97
                    parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmin
                    parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymax
                    parts(ig)%vx(1) = -parts(ig)%vx(1)
                    parts(ig)%vx(2) = -parts(ig)%vx(2)
                end if
                if ( (parts(i)%x(1)-vxmax)**2 + (parts(i)%x(2)-vymax)**2 < (scale_k*hsml)**2 .and. parts(i)%x(1) < vxmax) then
                    nghos_loc = nghos_loc + 1
                    ig = ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc
                    gind(nghos_loc) = i
                    parts(ig) = parts(i)
                    parts(ig)%indloc = ig
                    parts(ig)%itype = 97
                    parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmax
                    parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymax
                    parts(ig)%vx(1) = -parts(ig)%vx(1)
                    parts(ig)%vx(2) = -parts(ig)%vx(2)
                end if
                if ( (parts(i)%x(1)-vxmax)**2 + (parts(i)%x(2)-vymin)**2 < (scale_k*hsml)**2 .and. parts(i)%x(1) < vxmax) then
                    nghos_loc = nghos_loc + 1
                    ig = ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc
                    gind(nghos_loc) = i
                    parts(ig) = parts(i)
                    parts(ig)%indloc = ig
                    parts(ig)%itype = 97
                    parts(ig)%x(1) = -parts(ig)%x(1) + 2._f*vxmax
                    parts(ig)%x(2) = -parts(ig)%x(2) + 2._f*vymin
                    parts(ig)%vx(1) = -parts(ig)%vx(1)
                    parts(ig)%vx(2) = -parts(ig)%vx(2)
                end if
            end if
        end do
        
    end subroutine generate_ghost_part
    
    !==============================================================================================================================
    subroutine update_ghost_part
        
        implicit none
        integer:: i,j,k,d,ig,ir
        
        do i = 1,nghos_loc
            ig = ntotal_loc+nhalo_loc+nvirt_loc+i
            ir = gind(i)
            select case (parts(ig)%itype)
                case(99)
                    parts(ig)%rho = parts(ir)%rho
                    parts(ig)%p = parts(ir)%p
                    parts(ig)%vx(1) =-parts(ir)%vx(1)
                    parts(ig)%vx(2) = parts(ir)%vx(2)
                    parts(ig)%vx(3) = parts(ir)%vx(3)
                case(98)
                    parts(ig)%rho = parts(ir)%rho
                    parts(ig)%p = parts(ir)%p
                    parts(ig)%vx(1) = parts(ir)%vx(1)
                    parts(ig)%vx(2) =-parts(ir)%vx(2)
                    parts(ig)%vx(3) = parts(ir)%vx(3)
                case(97)
                    parts(ig)%rho = parts(ir)%rho
                    parts(ig)%p = parts(ir)%p
                    parts(ig)%vx(1) =-parts(ir)%vx(1)
                    parts(ig)%vx(2) =-parts(ir)%vx(2)
                    parts(ig)%vx(3) = parts(ir)%vx(3)
            end select
        end do
        
    end subroutine update_ghost_part
    
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
        if (ISNAN(beta)) beta = beta_max
        
        pv%rho = pr%rho
        pv%p = pr%p
        pv%vx(:) = (1._f-beta)*pr%vx(:)
        
    end subroutine virt_mirror

end module input_m
