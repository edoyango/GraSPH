module stress_update_m

    use datatypes, only: particles
    use param, only: dim,tenselem,f,pi,rh0,gamma,c,stress_mode
    
    integer,parameter:: Ide(tenselem) = [1,1,1,0,0,0], D2(tenselem) = [1,1,1,2,2,2]
    
contains
    
    !===============================================================================================================================
    subroutine stress_update(ki,dstrain,sig_0,part)
    
        implicit none
        integer,intent(in):: ki
        real(f),intent(in):: dstrain(tenselem),sig_0(tenselem)
        type(particles),intent(inout):: part
        real(f):: dpstrain(tenselem)=0._f
        
        select case (stress_mode)
        
            case ('SDP')
            
                call Drucker_Prager(dstrain,sig_0,part%sig,dpstrain)
                
            case ('EOS')
                
                call Equation_of_State(part%rho,part%sig)
                
            case default
            
                write(*,*) 'No stress mode selected!'
                stop
                
        end select
        
        part%p = -sum(part%sig(1:3))/3._f
        if (ki==4) part%pstrain(:) = part%pstrain(:) + dpstrain(:)
        
    end subroutine stress_update

    !===============================================================================================================================
    subroutine Equation_of_State(rho,sig)
    
        implicit none
        real(f),intent(in):: rho
        real(f),intent(out):: sig(tenselem)
        
        sig(1:3) = -rh0*c**2*((rho/rh0)**gamma-1._f)/gamma
        sig(4:tenselem) = 0._f
        
    end subroutine Equation_of_State
    
    !===============================================================================================================================
    subroutine Drucker_Prager(dstrain,sig_0,sig,dpstrain)
    
        implicit none
        real(f),intent(in):: dstrain(tenselem),sig_0(tenselem)
        real(f),intent(out):: sig(tenselem),dpstrain(tenselem)
        real(f),parameter:: phi=30._f,psi=0._f,coh=0._f,E=500000._f,v=0.3_f
        real(f),parameter:: D0 = E/((1._f+v)*(1._f-2._f*v)),Kb=E/(3._f*(1._f-2._f*v)),Gs=E/(2._f*(1._f+v))
        real(f),parameter:: DE(6,6) = D0*reshape([ 1._f-v,v     ,v     ,0._f       ,0._f       ,0._f,&
                                                   v     ,1._f-v,v     ,0._f       ,0._f       ,0._f,&
                                                   v     ,v     ,1._f-v,0._f       ,0._f       ,0._f,&
                                                   0._f  ,0._f  ,0._f  ,1._f-2._f*v,0._f       ,0._f,&
                                                   0._f  ,0._f  ,0._f  ,0._f       ,1._f-2._f*v,0._f,&
                                                   0._f  ,0._f  ,0._f  ,0._f       ,0._f       ,1._f-2._f*v ],[6,6])
        real(f),parameter:: sinphi = sin(phi*pi/180._f),sinpsi = tan(psi*pi/180._f),cosphi = cos(phi*pi/180._f)
        real(f),parameter:: alpha_phi = 2._f*sinphi/(sqrt(3._f)*(3._f-sinphi)),alpha_psi = 2._f*sinpsi/(sqrt(3._f)*(3._f-sinpsi)),&
            k_c = 6._f*coh*cosphi/(sqrt(3._f)*(3._f - sinphi))
        real(f):: dsig(tenselem),s(tenselem),dfdsig(tenselem),dgdsig(tenselem),I1,J2,dlambda,fy
        
        ! trial stress
        sig(:) = sig_0(:) + MATMUL(DE(:,:),dstrain(:))
        
        ! invariants and deviatoric stress tensor
        I1 = SUM(sig(1:dim))
        s(:) = sig(:) - I1/3._f*Ide(:)
        J2 = 0.5_f*SUM(s(:)*D2(:)*s(:))
        
        ! yield criterion
        fy = alpha_phi*I1 + SQRT(J2) - k_c
        
        if (fy > 0._f) then
            
            dfdsig(:) = alpha_phi*Ide(:) + s(:)/(2._f*sqrt(J2))
            dgdsig(:) = alpha_psi*Ide(:) + s(:)/(2._f*sqrt(J2))
            
            dlambda = fy/SUM(dfdsig(:)*D2*MATMUL(DE,dgdsig))
            
            dpstrain(:) = dlambda*dgdsig(:)
            sig(:) = sig(:) - MATMUL(DE,dpstrain(:))
            
            ! returning stress state to apex of yield surface
            I1 = SUM(sig(1:dim))
            if (I1 > k_c/alpha_phi) then
                sig(:) = k_c/alpha_phi/3._f*Ide(:)
            end if
            
        else
            
            sig(:) = sig(:)
            dpstrain(:) = 0._f
            
        end if
        
    end subroutine Drucker_Prager
    
end module stress_update_m
            
