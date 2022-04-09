module output_m
	
	use globvar, only: ntotal,nvirt,nghos,parts
	use param, only: output_directory,f,dim,tenselem
    
    use hdf5
    use h5lt
	
	integer,private:: i,d

contains

	!==============================================================================================================================
	subroutine output( ) 
	! Subroutine to write data to disk. Called intermittently, determined by save_step supplied at run-time
	
		use globvar,	only: itimestep,save_step
		use param,		only: output_phys,output_virt
		
		implicit none
        integer(HID_T):: file_id,real_group_id,virt_group_id,ghos_group_id,dspace_id,dset_id
        integer(HSIZE_T):: data_dims(2)
		integer:: n,ierr
		character(len=4)::number
        integer,allocatable:: idatatmp(:,:)
        real(f),allocatable:: fdatatmp(:,:)
		
		n = itimestep/save_step
		write(number,'(I4.4)') n

        ! Initializing hdf5 interface
        call h5open_f(ierr)
        
        ! Creating output file
        call h5fcreate_f(trim(output_directory)//"/sph_out"//number//".h5",h5F_ACC_TRUNC_F,file_id,ierr)
        
        ! Creating groups for each of real, virtual, and ghost particles
        call h5gcreate_f(file_id,"real",real_group_id,ierr)
        call h5gclose_f(real_group_id,ierr)
        
        call h5gcreate_f(file_id,"virt",virt_group_id,ierr)
        call h5gclose_f(virt_group_id,ierr)
        
        call h5gcreate_f(file_id,"ghos",ghos_group_id,ierr)
        call h5gclose_f(ghos_group_id,ierr)
        
        ! Writing particle data to real particles ----------------------------------------------------------------------------------
        ! Int data
        data_dims(:) = [ntotal,1]
        allocate(idatatmp(data_dims(1),data_dims(2)))
        
        ! particle index
        idatatmp(:,1) = parts(1:ntotal)%ind
        call H5LTmake_dataset_int_f(file_id,'real/ind',1,data_dims,idatatmp,ierr)
        
        ! process ID
        idatatmp(:,1) = 0
        call H5LTmake_dataset_int_f(file_id,'real/procid',1,data_dims,idatatmp,ierr)
        
        ! particle type
        idatatmp(:,1) = parts(1:ntotal)%itype
        call H5LTmake_dataset_int_f(file_id,'real/type',1,data_dims,idatatmp,ierr)
        deallocate(idatatmp)
        
        ! float data
        data_dims(:) = [dim,ntotal]
        allocate(fdatatmp(dim,ntotal))
        
        ! position
        do i = 1,ntotal
            fdatatmp(:,i) = parts(i)%x(:)
        end do
        call H5LTmake_dataset_double_f(file_id,'real/x',2,data_dims,fdatatmp,ierr)
        
        ! velocity
        do i = 1,ntotal
            fdatatmp(:,i) = parts(i)%vx(:)
        end do
        call H5LTmake_dataset_double_f(file_id,'real/v',2,data_dims,fdatatmp,ierr)
        
        deallocate(fdatatmp)
        
        data_dims(:) = [tenselem,ntotal]
        allocate(fdatatmp(tenselem,ntotal))
        do i = 1,ntotal
            fdatatmp(:,i) = parts(i)%strain(:)
        end do
        call H5LTmake_dataset_double_f(file_id,'real/strain',2,data_dims,fdatatmp,ierr)
        do i = 1,ntotal
            fdatatmp(:,i) = parts(i)%sig(:)
        end do
        call H5LTmake_dataset_double_f(file_id,'real/stress',2,data_dims,fdatatmp,ierr)
        do i = 1,ntotal
            fdatatmp(:,i) = parts(i)%pstrain(:)
        end do
        call H5LTmake_dataset_double_f(file_id,'real/pstrain',2,data_dims,fdatatmp,ierr)
        
        deallocate(fdatatmp)
        
        ! density
        data_dims(:) = [ntotal,1]
        allocate(fdatatmp(ntotal,1))
        fdatatmp(:,1) = parts(1:ntotal)%rho
        call H5LTmake_dataset_double_f(file_id,'real/rho',1,data_dims,fdatatmp,ierr)
        
        ! pressure
        fdatatmp(:,1) = parts(1:ntotal)%p
        call H5LTmake_dataset_double_f(file_id,'real/p',1,data_dims,fdatatmp,ierr)
        
        deallocate(fdatatmp)
        
        ! Writing particle data to virtual particles -------------------------------------------------------------------------------
        ! Int data
        data_dims(:) = [nvirt,1]
        allocate(idatatmp(data_dims(1),data_dims(2)))
        
        ! particle index
        idatatmp(:,1) = parts(ntotal+1:ntotal+nvirt)%ind
        call H5LTmake_dataset_int_f(file_id,'virt/ind',1,data_dims,idatatmp,ierr)
        
        ! process ID
        idatatmp(:,1) = 0
        call H5LTmake_dataset_int_f(file_id,'virt/procid',1,data_dims,idatatmp,ierr)
        
        ! particle type
        idatatmp(:,1) = parts(ntotal+1:ntotal+nvirt)%itype
        call H5LTmake_dataset_int_f(file_id,'virt/type',1,data_dims,idatatmp,ierr)
        deallocate(idatatmp)
        
        ! float data
        data_dims(:) = [dim,nvirt]
        allocate(fdatatmp(dim,nvirt))
        
        ! position
        do i = 1,nvirt
            fdatatmp(:,i) = parts(ntotal+i)%x(:)
        end do
        call H5LTmake_dataset_double_f(file_id,'virt/x',2,data_dims,fdatatmp,ierr)
        
        ! velocity
        do i = 1,nvirt
            fdatatmp(:,i) = parts(ntotal+i)%vx(:)
        end do
        call H5LTmake_dataset_double_f(file_id,'virt/v',2,data_dims,fdatatmp,ierr)
        
        deallocate(fdatatmp)
        
        ! density
        data_dims(:) = [nvirt,1]
        allocate(fdatatmp(nvirt,1))
        fdatatmp(:,1) = parts(ntotal+1:ntotal+nvirt)%rho
        call H5LTmake_dataset_double_f(file_id,'virt/rho',1,data_dims,fdatatmp,ierr)
        
        ! pressure
        fdatatmp(:,1) = parts(ntotal+1:ntotal+nvirt)%p
        call H5LTmake_dataset_double_f(file_id,'virt/p',1,data_dims,fdatatmp,ierr)
        
        deallocate(fdatatmp)
        
        ! Writing particle data to ghost particles ---------------------------------------------------------------------------------
        ! Int data
        data_dims(:) = [nghos,1]
        allocate(idatatmp(data_dims(1),data_dims(2)))
        
        ! particle index
        idatatmp(:,1) = parts(ntotal+nvirt+1:ntotal+nvirt+nghos)%ind
        call H5LTmake_dataset_int_f(file_id,'ghos/ind',1,data_dims,idatatmp,ierr)
        
        ! process ID
        idatatmp(:,1) = 0
        call H5LTmake_dataset_int_f(file_id,'ghos/procid',1,data_dims,idatatmp,ierr)
        
        ! particle type
        idatatmp(:,1) = parts(ntotal+nvirt+1:ntotal+nvirt+nghos)%itype
        call H5LTmake_dataset_int_f(file_id,'ghos/type',1,data_dims,idatatmp,ierr)
        deallocate(idatatmp)
        
        ! float data
        data_dims(:) = [dim,nghos]
        allocate(fdatatmp(dim,nghos))
        
        ! position
        do i = 1,nghos
            fdatatmp(:,i) = parts(ntotal+nvirt+i)%x(:)
        end do
        call H5LTmake_dataset_double_f(file_id,'ghos/x',2,data_dims,fdatatmp,ierr)
        
        ! velocity
        do i = 1,nghos
            fdatatmp(:,i) = parts(ntotal+nvirt+i)%vx(:)
        end do
        call H5LTmake_dataset_double_f(file_id,'ghos/v',2,data_dims,fdatatmp,ierr)
        
        deallocate(fdatatmp)
        
        ! density
        data_dims(:) = [nghos,1]
        allocate(fdatatmp(nghos,1))
        fdatatmp(:,1) = parts(ntotal+nvirt+1:ntotal+nvirt+nghos)%rho
        call H5LTmake_dataset_double_f(file_id,'ghos/rho',1,data_dims,fdatatmp,ierr)
        
        ! pressure
        fdatatmp(:,1) = parts(ntotal+nvirt+1:ntotal+nvirt+nghos)%p
        call H5LTmake_dataset_double_f(file_id,'ghos/p',1,data_dims,fdatatmp,ierr)
        
        deallocate(fdatatmp)
        
        call h5fclose_f(file_id,ierr)
        
	end subroutine output
	
	!==============================================================================================================================
	subroutine write_ini_config
	! Subroutine for writing initial configuration data using intrinsic IO
	
		implicit none
        integer(HID_T):: file_id,real_group_id,virt_group_id,ghos_group_id,dspace_id,dset_id
        integer(HSIZE_T):: data_dims(2)
		integer:: ierr
        integer,allocatable:: idatatmp(:,:)
        real(f),allocatable:: fdatatmp(:,:)
        
        ! Initializing hdf5 interface
        call h5open_f(ierr)
        
        ! Creating output file
        call h5fcreate_f(trim(output_directory)//"/sph_out0000.h5",h5F_ACC_TRUNC_F,file_id,ierr)
        
        ! Creating groups for each of real, virtual, and ghost particles
        call h5gcreate_f(file_id,"real",real_group_id,ierr)
        call h5gclose_f(real_group_id,ierr)
        
        call h5gcreate_f(file_id,"virt",virt_group_id,ierr)
        call h5gclose_f(virt_group_id,ierr)
        
        call h5gcreate_f(file_id,"ghos",ghos_group_id,ierr)
        call h5gclose_f(ghos_group_id,ierr)
        
        ! Writing particle data to real particles ----------------------------------------------------------------------------------
        ! Int data
        data_dims(:) = [ntotal,1]
        allocate(idatatmp(data_dims(1),data_dims(2)))
        
        ! particle index
        idatatmp(:,1) = parts(1:ntotal)%ind
        call H5LTmake_dataset_int_f(file_id,'real/ind',1,data_dims,idatatmp,ierr)
        
        ! process ID
        idatatmp(:,1) = 0
        call H5LTmake_dataset_int_f(file_id,'real/procid',1,data_dims,idatatmp,ierr)
        
        ! particle type
        idatatmp(:,1) = parts(1:ntotal)%itype
        call H5LTmake_dataset_int_f(file_id,'real/type',1,data_dims,idatatmp,ierr)
        deallocate(idatatmp)
        
        ! float data
        data_dims(:) = [dim,ntotal]
        allocate(fdatatmp(dim,ntotal))
        
        ! position
        do i = 1,ntotal
            fdatatmp(:,i) = parts(i)%x(:)
        end do
        call H5LTmake_dataset_double_f(file_id,'real/x',2,data_dims,fdatatmp,ierr)
        
        ! velocity
        do i = 1,ntotal
            fdatatmp(:,i) = parts(i)%vx(:)
        end do
        call H5LTmake_dataset_double_f(file_id,'real/v',2,data_dims,fdatatmp,ierr)
        
        deallocate(fdatatmp)
        
        data_dims(:) = [tenselem,ntotal]
        allocate(fdatatmp(tenselem,ntotal))
        do i = 1,ntotal
            fdatatmp(:,i) = parts(i)%strain(:)
        end do
        call H5LTmake_dataset_double_f(file_id,'real/strain',2,data_dims,fdatatmp,ierr)
        do i = 1,ntotal
            fdatatmp(:,i) = parts(i)%sig(:)
        end do
        call H5LTmake_dataset_double_f(file_id,'real/stress',2,data_dims,fdatatmp,ierr)
        
        deallocate(fdatatmp)
        
        ! density
        data_dims(:) = [ntotal,1]
        allocate(fdatatmp(ntotal,1))
        fdatatmp(:,1) = parts(1:ntotal)%rho
        call H5LTmake_dataset_double_f(file_id,'real/rho',1,data_dims,fdatatmp,ierr)
        
        ! pressure
        fdatatmp(:,1) = parts(1:ntotal)%p
        call H5LTmake_dataset_double_f(file_id,'real/p',1,data_dims,fdatatmp,ierr)
        
        deallocate(fdatatmp)
		
	end subroutine write_ini_config
    
end module output_m
