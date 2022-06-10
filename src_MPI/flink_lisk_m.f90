module flink_list_m
    
    use datatypes,        only: particles
    use globvar,        only: ntotal_loc,nhalo_loc,nvirt_loc,nghos_loc,parts,scale_k,niac,pairs,maxinter
    use globvar_para,    only: procid,numprocs
    use param,            only: dim,hsml,f
    
    use error_msg_m,    only: error_msg
    use kernel_m,        only: kernel
    
    type particleincellarray
        type(particles),pointer:: p
    end type particleincellarray
    
    integer,private:: i,j,k,d,icell,jcell,kcell,xi,yi,zi
    real(f),private:: dcell
    
    public:: flink_list
    private:: particleincellarray,check_if_interact
        
contains
    
    !==============================================================================================================================
    subroutine flink_list( )
    ! save as above, but for 3D
        
        implicit none
        integer,parameter:: maxpcell=125
        integer:: ngridx(3)
        real(f):: mingridx(3),maxgridx(3)
        integer,allocatable:: pincell(:,:,:)
        type(particleincellarray),allocatable:: cells(:,:,:,:)
        integer,parameter:: sweep(3,13) = reshape((/-1,-1,-1,&
                                                    -1,-1, 0,&
                                                    -1,-1, 1,&
                                                    -1, 0,-1,&
                                                    -1, 0, 0,&
                                                    -1, 0, 1,&
                                                    -1, 1,-1,&
                                                    -1, 1, 0,&
                                                    -1, 1, 1,&
                                                     0,-1,-1,&
                                                     0,-1, 0,&
                                                     0,-1, 1,&
                                                     0, 0,-1/),(/3,13/))
        
        !Determining bounding box extents
        mingridx(:) = parts(1)%x(:)
        maxgridx(:) = parts(1)%x(:)
        do i = 2,ntotal_loc
            do d = 1,dim
                mingridx(d) = MIN(mingridx(d),parts(i)%x(d))
                maxgridx(d) = MAX(maxgridx(d),parts(i)%x(d))
            end do
        end do
        
        !Determining number of grid cells in each direction
        dcell = scale_k*hsml
        maxgridx(:) = maxgridx(:) + 3._f*dcell
        mingridx(:) = mingridx(:) - 3._f*dcell
        ngridx(:) = int((maxgridx(:) - mingridx(:))/dcell) + 1
        maxgridx(:) = mingridx(:) + ngridx(:)*dcell
        
        allocate( pincell(ngridx(1),ngridx(2),ngridx(3)),&
                    cells(maxpcell,ngridx(1),ngridx(2),ngridx(3)) )
                    
        pincell(:,:,:) = 0
        
        do i=1,ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc
            if (.not.any(parts(i)%x(:) < mingridx(:) .or. parts(i)%x(:) >= maxgridx(:))) then
                icell = int((parts(i)%x(1) - mingridx(1))/dcell) + 1
                jcell = int((parts(i)%x(2) - mingridx(2))/dcell) + 1 
                kcell = int((parts(i)%x(3) - mingridx(3))/dcell) + 1 
                pincell(icell,jcell,kcell) = pincell(icell,jcell,kcell) + 1
                cells(pincell(icell,jcell,kcell),icell,jcell,kcell)%p => parts(i)
            end if
        enddo
        
        niac = 0
        do icell = 2,ngridx(1)-2
            do jcell = 2,ngridx(2)-2
                do kcell = 2,ngridx(3)-2
                    if (pincell(icell,jcell,kcell) > 0) then
                    
                        ! finding pairs within cell icell,jcell
                        do i = 1,pincell(icell,jcell,kcell)-1
                            do j = i+1,pincell(icell,jcell,kcell)
                                call check_if_interact(cells(i,icell,jcell,kcell)%p,cells(j,icell,jcell,kcell)%p)
                            end do
                        end do
                        
                        ! finding pairs between particles in cell icell,jcell and particles in cell xi,yi
                        do k = 1,13
                            xi = icell + sweep(1,k)
                            yi = jcell + sweep(2,k)
                            zi = kcell + sweep(3,k)
                            do i = 1,pincell(icell,jcell,kcell)
                                do j = 1,pincell(xi,yi,zi)
                                    call check_if_interact(cells(i,icell,jcell,kcell)%p,cells(j,xi,yi,zi)%p)
                                end do
                            end do
                        end do
                    
                    endif
                end do
            end do
        end do

    end subroutine flink_list
    
    !==============================================================================================================================
    subroutine check_if_interact(p_i,p_j)
    ! subroutine to chekc if two particles are interacting and consequently adding to pair list
        
        implicit none
        type(particles),intent(in),target:: p_i,p_j
        real(f):: dxiac(dim),r
        
        ! only consider interactions when real-real are involved
        if ( (p_i%itype.eq.1 .or. p_j%itype.eq.1) ) then
            dxiac(:) = p_i%x(:) - p_j%x(:)
            r = SQRT(SUM(dxiac*dxiac))
            if (r <= hsml*scale_k) then
                niac = niac + 1
                if (niac < maxinter) then
                    pairs(niac)%i => p_i
                    pairs(niac)%j => p_j
                    call kernel(r,dxiac,hsml,pairs(niac)%w,pairs(niac)%dwdx(:))
                else
                    print *,' >>> Error <<< : Too many interactions'
                    stop
                end if
            end if
        end if
        
    end subroutine check_if_interact
    
end module flink_list_m
