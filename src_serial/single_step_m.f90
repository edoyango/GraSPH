module single_step_m

    public:: single_step

contains

    !==============================================================================================================================
    subroutine single_step(ki,dvxdti,drhoi) 
    ! Container subroutine for all the rate-of-change calculations. Rate-of-changes are calculated seperately and then summed as
    ! required
        
        use globvar, only: ntotal,nvirt,nghos,niac,pairs,parts
        use param, only: dim,f,g
        
        use input_m, only: virt_mirror
        use material_rates_m
        
        implicit none
        integer,intent(in):: ki
        real(f),intent(out):: dvxdti(dim,ntotal),drhoi(ntotal)
        integer:: i,j,k,d
        real(f):: t1,t2
        real(f),allocatable:: indvxdt(:,:),ardvxdt(:,:),exdvxdt(:,:),cdrhodt(:)
        
        allocate( indvxdt(dim,ntotal+nvirt+nghos),&
                ardvxdt(dim,ntotal+nvirt+nghos),&
                exdvxdt(dim,ntotal+nvirt+nghos),&
                cdrhodt(ntotal+nvirt+nghos) )
                
        cdrhodt(1:ntotal) = 0_f
        indvxdt(:,1:ntotal) = 0_f
        ardvxdt(:,1:ntotal) = 0_f
        exdvxdt(1:dim-1,1:ntotal) = 0_f
        exdvxdt(dim,1:ntotal) = -g
        
        ! looping through interaction pairs to calculate forces/density change
        do k = 1,niac
            i = pairs(k)%i
            j = pairs(k)%j
            if (parts(i)%itype > 0 .and. parts(j)%itype < 0) then
                call virt_mirror(parts(i),parts(j))
            elseif (parts(i)%itype < 0 .and. parts(j)%itype > 0) then
                call virt_mirror(parts(j),parts(i))
            end if
                
                !Density approximation or change rate
                call con_density(ki,parts(i),parts(j),pairs(k)%dwdx,cdrhodt)
                
                !Internal force due to pressure
                call int_force(ki,parts(i),parts(j),pairs(k)%dwdx,indvxdt)
                
                !Artificial viscosity:
                call art_visc(ki,parts(i),parts(j),pairs(k)%dwdx,ardvxdt)
                
        end do
        
        !Convert velocity, force, and energy to f and dfdt
        dvxdti(1:dim,1:ntotal) = indvxdt(1:dim,1:ntotal) + exdvxdt(1:dim,1:ntotal) + ardvxdt(1:dim,1:ntotal)
        drhoi(1:ntotal) = cdrhodt(1:ntotal)
        
        deallocate( indvxdt,ardvxdt,exdvxdt,cdrhodt )
    
    end subroutine single_step

end module single_step_m
