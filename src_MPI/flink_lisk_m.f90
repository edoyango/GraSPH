module flink_list_m
    
    use datatypes,        only: particles
    use globvar,        only: ntotal_loc,nhalo_loc,nvirt_loc,nghos_loc,parts,scale_k,niac,pairs,maxinter
    use globvar_para,    only: procid,numprocs
    use param,            only: dim,hsml,f
    
    use error_msg_m,    only: error_msg
    use kernel_m,        only: kernel
    
    public:: flink_list
    private:: check_if_interact
        
contains

    !==============================================================================================================================
    subroutine flink_list( )
    ! save as above, but for 3D
        
        implicit none
        integer,parameter:: maxpcell=125
        integer:: ngridx(3),jth,i,j,k,d,icell,jcell,kcell,xi,yi,zi
        real(f):: mingridx(3),maxgridx(3),dcell
        integer,allocatable:: pincell(:,:,:),gridind(:,:),cells(:,:,:,:)
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
        do i = 2,ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc
            do d = 1,dim
                mingridx(d) = MIN(mingridx(d),parts(i)%x(d))
                maxgridx(d) = MAX(maxgridx(d),parts(i)%x(d))
            end do
        end do
        
        !Determining number of grid cells in each direction
        dcell = scale_k*hsml
        maxgridx(:) = maxgridx(:) + 2._f*dcell
        mingridx(:) = mingridx(:) - 2._f*dcell
        ngridx(:) = int((maxgridx(:) - mingridx(:))/dcell) + 1
        maxgridx(:) = mingridx(:) + ngridx(:)*dcell
        
        allocate( pincell(ngridx(1),ngridx(2),ngridx(3)),&
                    gridind(dim,ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc),&
                    cells(maxpcell,ngridx(1),ngridx(2),ngridx(3)) )
                    
        pincell(:,:,:) = 0
        
        do i=1,ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc
            gridind(:,i) = int((parts(i)%x(:) - mingridx(:))/dcell) + 1
            icell = gridind(1,i)
            jcell = gridind(2,i) 
            kcell = gridind(3,i) 
            pincell(icell,jcell,kcell) = pincell(icell,jcell,kcell) + 1
            cells(pincell(icell,jcell,kcell),icell,jcell,kcell) = i
        enddo
        
        niac = 0
        do i = 1,ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc
            icell = gridind(1,i); jcell = gridind(2,i); kcell = gridind(3,i)
            do j = 1,pincell(icell,jcell,kcell)
                jth = cells(j,icell,jcell,kcell)
                if (jth > i) then
                    call check_if_interact(parts(i),parts(jth))
                end if
            end do
            do k = 1,13
                xi = icell + sweep(1,k)
                yi = jcell + sweep(2,k)
                zi = kcell + sweep(3,k)
                do j = 1,pincell(xi,yi,zi)
                    jth = cells(j,xi,yi,zi)
                    call check_if_interact(parts(i),parts(jth))
                end do
            end do
        end do

    end subroutine flink_list
    
    !==============================================================================================================================
    subroutine check_if_interact(p_i,p_j)
    ! subroutine to chekc if two particles are interacting and consequently adding to pair list
        
        implicit none
        type(particles),intent(in):: p_i,p_j
        real(f):: dxiac(dim),r
        
        ! only consider interactions when real-real are involved
        if ( (p_i%itype.eq.1 .or. p_j%itype.eq.1) ) then
            dxiac(:) = p_i%x(:) - p_j%x(:)
            r = SQRT(SUM(dxiac*dxiac))
            if (r <= hsml*scale_k) then
                niac = niac + 1
                if (niac < maxinter) then
                    pairs(niac)%i = p_i%indloc
                    pairs(niac)%j = p_j%indloc
                    call kernel(r,dxiac,hsml,pairs(niac)%w,pairs(niac)%dwdx(:))
                else
                    print *,' >>> Error <<< : Too many interactions'
                    stop
                end if
            end if
        end if
        
    end subroutine check_if_interact
    
end module flink_list_m
