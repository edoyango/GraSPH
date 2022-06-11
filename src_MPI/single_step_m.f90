module single_step_m

    use input_m, only: update_ghost_part,virt_mirror
    
    public:: single_step
    
contains

    !==============================================================================================================================
    subroutine single_step(ki,dvxdti,drhoi)
    ! Container subroutine for all the rate-of-change calculations. Rate-of-changes are calculated seperately and then summed as
    ! required
        
        use globvar,             only: ntotal_loc,nhalo_loc,nvirt_loc,nghos_loc,parts,pairs,t_dist,niac
        use mpi_f08,                 only: MPI_WTIME
        use param,                 only: dim,rh0,c,gamma,f,g
        
        use material_rates_m,    only: int_force,art_visc,con_density,ext_force
        use ORB_sr_m,             only: ORB_sendrecv_haloupdate
        
        implicit none
        integer,intent(in):: ki
        real(f),intent(out):: dvxdti(dim,ntotal_loc),drhoi(ntotal_loc)
        integer:: i,j,k
        real(f),allocatable:: indvxdt(:,:),ardvxdt(:,:),exdvxdt(:,:),codrhodt(:)
        
        t_dist = t_dist - MPI_WTIME()
        if (ki.ne.1) call ORB_sendrecv_haloupdate(ki)
        t_dist = t_dist + MPI_WTIME()
        
        parts(1:ntotal_loc+nhalo_loc)%p = rh0*c**2*((parts(1:ntotal_loc+nhalo_loc)%rho/rh0)**gamma-1_f)/gamma
        
        if (ki.ne.1) call update_ghost_part
        
        allocate(indvxdt(dim,ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc),&
            ardvxdt(dim,ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc),&
            exdvxdt(dim,ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc),&
            codrhodt(ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc) )
        
        codrhodt(1:ntotal_loc) = 0._f
        indvxdt(:,1:ntotal_loc) = 0._f
        ardvxdt(:,1:ntotal_loc) = 0._f
        exdvxdt(1:dim-1,1:ntotal_loc) = 0._f
        exdvxdt(dim,1:ntotal_loc) = -g
        
        do k = 1,niac
            
            i = pairs(k)%i
            j = pairs(k)%j
            
            if (parts(i)%itype > 0 .and. parts(j)%itype < 0) then
                call virt_mirror(parts(i),parts(j))
            elseif (parts(i)%itype < 0 .and. parts(j)%itype > 0) then
                call virt_mirror(parts(j),parts(i))
            end if
            
            !Internal force due to pressure
            call int_force(ki,parts(i),parts(j),pairs(k)%dwdx,indvxdt(:,i),indvxdt(:,j))
            
            !Artificial viscosity:
            call art_visc(ki,parts(i),parts(j),pairs(k)%dwdx,ardvxdt(:,i),ardvxdt(:,j))
            
            !Density approximation or change rate
            call con_density(ki,parts(i),parts(j),pairs(k)%dwdx,codrhodt(i),codrhodt(j))
            
        end do
        
        !Convert velocity, force, and energy to f and dfdt
        dvxdti(1:dim,1:ntotal_loc) = indvxdt(1:dim,1:ntotal_loc) + exdvxdt(1:dim,1:ntotal_loc) + ardvxdt(1:dim,1:ntotal_loc)
        drhoi(1:ntotal_loc) = codrhodt(1:ntotal_loc)
        
        deallocate( indvxdt,ardvxdt,exdvxdt,codrhodt )
        
    end subroutine single_step

end module single_step_m
