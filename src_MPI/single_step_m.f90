module single_step_m

    use input_m, only: update_ghost_part,virt_mirror
    
    public:: single_step
    
contains

    !==============================================================================================================================
    subroutine single_step(ki,dvxdti,drhoi)
    ! Container subroutine for all the rate-of-change calculations. Rate-of-changes are calculated seperately and then summed as
    ! required
        
        use globvar,             only: ntotal_loc,nhalo_loc,nvirt_loc,nghos_loc,parts,pairs,t_dist,niac
        use mpi,                 only: MPI_WTIME
        use param,                 only: dim,rh0,c,gamma,f,g
        
        use material_rates_m,    only: int_force,art_visc,con_density,ext_force
        use ORB_sr_m,             only: ORB_sendrecv_haloupdate
        
        implicit none
        integer,intent(in):: ki
        real(f),intent(out):: dvxdti(dim,ntotal_loc),drhoi(ntotal_loc)
        integer:: k
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
                
            if (pairs(k)%i%itype > 0 .and. pairs(k)%j%itype < 0) then
                call virt_mirror(pairs(k)%i,pairs(k)%j)
            elseif (pairs(k)%i%itype < 0 .and. pairs(k)%j%itype > 0) then
                call virt_mirror(pairs(k)%j,pairs(k)%i)
            end if
            
            !Internal force due to pressure
            call int_force(ki,pairs(k),indvxdt)
            
            !Artificial viscosity:
            call art_visc(ki,pairs(k),ardvxdt)
            
            !Density approximation or change rate
            call con_density(ki,pairs(k),codrhodt)
            
        end do
        
        !Convert velocity, force, and energy to f and dfdt
        dvxdti(1:dim,1:ntotal_loc) = indvxdt(1:dim,1:ntotal_loc) + exdvxdt(1:dim,1:ntotal_loc) + ardvxdt(1:dim,1:ntotal_loc)
        drhoi(1:ntotal_loc) = codrhodt(1:ntotal_loc)
        
        deallocate( indvxdt,ardvxdt,exdvxdt,codrhodt )
        
    end subroutine single_step

end module single_step_m
