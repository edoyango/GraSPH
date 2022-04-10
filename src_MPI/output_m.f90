module output_m

	use globvar,		only: ntotal_loc,nhalo_loc,nvirt_loc,nghos_loc,parts,itimestep,save_step,ntotal,nvirt
	use globvar_para,	only: procid,numprocs,ierr,MPI_ftype
	use param,			only: f,dim,output_directory,output_phys,output_halo,output_virt,output_flt_type,tenselem
    
    use hdf5
    use mpi
    
    use hdf5_parallel_io_helper_m, only: hdf5_parallel_fileopen,hdf5_parallel_write
	
	implicit none
	character(len=220),private:: filepath
	character(len=4),private::number
	integer,allocatable,private:: int_tmp(:,:)
	real(f),allocatable,private:: dbl_tmp(:,:)
	integer,private:: n,i,j,k,d,request(4),status(MPI_STATUS_SIZE,4),pos0,pos1
    integer(HID_T),private:: fid,real_group_id,virt_group_id,ghos_group_id
    integer(HSIZE_T),private:: global_dims(2)
    integer(HSSIZE_T),private:: displ(2)
	
	public:: output,write_ini_config
	
contains

	!==============================================================================================================================
	subroutine output
	! Subroutine to write data to disk. Called intermittently, determined by save_step supplied at run-time
		
		implicit none
        integer:: ntotal_glob(numprocs),nhalo_glob(numprocs),nvirt_glob(numprocs),nghos_glob(numprocs)
        
		! Exchanging how many particles each process will output data for
		call MPI_IALLGATHER(ntotal_loc,1,MPI_INTEGER,ntotal_glob,1,MPI_INTEGER,MPI_COMM_WORLD,request(1),ierr)
		call MPI_IALLGATHER(nhalo_loc,1,MPI_INTEGER,nhalo_glob,1,MPI_INTEGER,MPI_COMM_WORLD,request(2),ierr)
		call MPI_IALLGATHER(nvirt_loc,1,MPI_INTEGER,nvirt_glob,1,MPI_INTEGER,MPI_COMM_WORLD,request(3),ierr)
		call MPI_IALLGATHER(nghos_loc,1,MPI_INTEGER,nghos_glob,1,MPI_INTEGER,MPI_COMM_WORLD,request(4),ierr)
		
		!Format number
		n = itimestep/save_step
		write(number,'(I4.4)') n
        
        ! initializing hdf5
        call h5open_f(ierr)
        
        ! creating hdf5 output file (hdf5_parallel_fileopen is a custom hdf5 wrapper subroutine)
        filepath = trim(output_directory)//"/sph_out"//number//".h5"
        call hdf5_parallel_fileopen(filepath,fid)
        
        ! creating groups in hdf5 output for real, halo, virtual, and ghost particles
        call h5gcreate_f(fid,"real",real_group_id,ierr)
        call h5gclose_f(real_group_id,ierr)
        
        call h5gcreate_f(fid,"halo",ghos_group_id,ierr)
        call h5gclose_f(ghos_group_id,ierr)
        
        call h5gcreate_f(fid,"virt",virt_group_id,ierr)
        call h5gclose_f(virt_group_id,ierr)
        
        call h5gcreate_f(fid,"ghos",ghos_group_id,ierr)
        call h5gclose_f(ghos_group_id,ierr)
		
        ! Waiting for non-blocking exchange to complete
		call MPI_WAITALL(4,request,status,ierr)
        
        ! Writing data for real particles ------------------------------------------------------------------------------------------
        ! defining array shapes and displacments
        global_dims(:) = [dim,ntotal]
        displ(:) = [0,SUM(ntotal_glob(1:procid))] ! hdf5 works with 0-indexing
        
        allocate( int_tmp(ntotal_loc,1),dbl_tmp(ntotal_loc,1) )
        
        ! writing 1D arrays
        pos0 = 1
        pos1 = ntotal_loc
        int_tmp(:,1) = parts(pos0:pos1)%indglob
        call hdf5_parallel_write(fid,'real/ind',displ(2),global_dims(2),int_tmp(:,1))
        int_tmp(:,1) = procid
        call hdf5_parallel_write(fid,'real/procid',displ(2),global_dims(2),int_tmp(:,1))
        int_tmp(:,1) = parts(pos0:pos1)%itype
        call hdf5_parallel_write(fid,'real/type',displ(2),global_dims(2),int_tmp(:,1))
        dbl_tmp(:,1) = parts(pos0:pos1)%rho
        call hdf5_parallel_write(fid,'real/rho',displ(2),global_dims(2),dbl_tmp(:,1))
        dbl_tmp(:,1) = parts(pos0:pos1)%p
        call hdf5_parallel_write(fid,'real/p',displ(2),global_dims(2),dbl_tmp(:,1))
        
        deallocate( int_tmp,dbl_tmp )
        
        ! Converting position data to contiguous array
        allocate(dbl_tmp(dim,ntotal_loc))
        do i = 1,ntotal_loc
            dbl_tmp(:,i) = parts(i)%x(:)
        end do
        
        ! writing position
        call hdf5_parallel_write(fid,'real/x',displ,global_dims,dbl_tmp)
        
        ! Converting velocity data to contiguous array
        do i = 1,ntotal_loc
            dbl_tmp(:,i) = parts(i)%vx(:)
        end do
        
        ! writing velocity
        call hdf5_parallel_write(fid,'real/v',displ,global_dims,dbl_tmp)
        
        ! cleanup
        deallocate(dbl_tmp)
        
        ! Writing tensor data
        ! Converting data to contiguous array
        allocate(dbl_tmp(tenselem,ntotal_loc))
        
        global_dims(:) = [tenselem,ntotal]
        
        do i = 1,ntotal_loc
            dbl_tmp(:,i) = parts(i)%strain(:)
        end do
        call hdf5_parallel_write(fid,'real/strain',displ,global_dims,dbl_tmp)
        do i = 1,ntotal_loc
            dbl_tmp(:,i) = parts(i)%pstrain(:)
        end do
        call hdf5_parallel_write(fid,'real/pstrain',displ,global_dims,dbl_tmp)
        do i = 1,ntotal_loc
            dbl_tmp(:,i) = parts(i)%sig(:)
        end do
        call hdf5_parallel_write(fid,'real/stress',displ,global_dims,dbl_tmp)
        
        deallocate(dbl_tmp)
        
        ! Writing data for halo particles ------------------------------------------------------------------------------------------
        ! defining array shapes and displacments
        global_dims(:) = [dim,sum(nhalo_glob)]
        displ(:) = [0,SUM(nhalo_glob(1:procid))] ! hdf5 works with 0-indexing
        
        allocate( int_tmp(nhalo_loc,1),dbl_tmp(nhalo_loc,1) )
        
        ! writing 1D arrays
        pos0 = ntotal_loc+1
        pos1 = ntotal_loc+nhalo_loc
        int_tmp(:,1) = parts(pos0:pos1)%indglob
        call hdf5_parallel_write(fid,'halo/ind',displ(2),global_dims(2),int_tmp(:,1))
        int_tmp(:,1) = procid
        call hdf5_parallel_write(fid,'halo/procid',displ(2),global_dims(2),int_tmp(:,1))
        int_tmp(:,1) = parts(pos0:pos1)%itype
        call hdf5_parallel_write(fid,'halo/type',displ(2),global_dims(2),int_tmp(:,1))
        dbl_tmp(:,1) = parts(pos0:pos1)%rho
        call hdf5_parallel_write(fid,'halo/rho',displ(2),global_dims(2),dbl_tmp(:,1))
        dbl_tmp(:,1) = parts(pos0:pos1)%p
        call hdf5_parallel_write(fid,'halo/p',displ(2),global_dims(2),dbl_tmp(:,1))
        
        deallocate( int_tmp,dbl_tmp )
        
        ! Converting position data to contiguous array
        allocate(dbl_tmp(dim,nhalo_loc))
        do i = 1,nhalo_loc
            dbl_tmp(:,i) = parts(ntotal_loc+i)%x(:)
        end do
        
        ! writing position
        call hdf5_parallel_write(fid,'halo/x',displ,global_dims,dbl_tmp)
        
        ! Converting velocity data to contiguous array
        do i = 1,nhalo_loc
            dbl_tmp(:,i) = parts(ntotal_loc+i)%vx(:)
        end do

        ! writing velocity
        call hdf5_parallel_write(fid,'halo/v',displ,global_dims,dbl_tmp)
        
        ! cleanup
        deallocate(dbl_tmp)
        
        allocate(dbl_tmp(tenselem,nhalo_loc))
        
        global_dims(:) = [tenselem,sum(nhalo_glob)]
        
        do i = 1,nhalo_loc
            dbl_tmp(:,i) = parts(ntotal_loc+i)%sig(:)
        end do
        call hdf5_parallel_write(fid,'halo/stress',displ,global_dims,dbl_tmp)
        
        deallocate(dbl_tmp)
        
        ! Writing data for virtual particles ---------------------------------------------------------------------------------------
        ! defining array shapes and displacments
        global_dims(:) = [dim,sum(nvirt_glob)]
        displ(:) = [0,SUM(nvirt_glob(1:procid))] ! hdf5 works with 0-indexing
        
        allocate(int_tmp(nvirt_loc,1),dbl_tmp(nvirt_loc,1))
        
        ! writing 1D arrays
        pos0 = ntotal_loc+nhalo_loc+1
        pos1 = ntotal_loc+nhalo_loc+nvirt_loc
        int_tmp(:,1) = parts(pos0:pos1)%indglob
        call hdf5_parallel_write(fid,'virt/ind',displ(2),global_dims(2),int_tmp(:,1))
        int_tmp(:,1) = procid
        call hdf5_parallel_write(fid,'virt/procid',displ(2),global_dims(2),int_tmp(:,1))
        int_tmp(:,1) = parts(pos0:pos1)%itype
        call hdf5_parallel_write(fid,'virt/type',displ(2),global_dims(2),int_tmp(:,1))
        dbl_tmp(:,1) = parts(pos0:pos1)%rho
        call hdf5_parallel_write(fid,'virt/rho',displ(2),global_dims(2),dbl_tmp(:,1))
        dbl_tmp(:,1) = parts(pos0:pos1)%p
        call hdf5_parallel_write(fid,'virt/p',displ(2),global_dims(2),dbl_tmp(:,1))
        
        deallocate( int_tmp,dbl_tmp)
        
        ! Converting position data to contiguous array
        allocate(dbl_tmp(dim,nvirt_loc))
        do i = 1,nvirt_loc
            dbl_tmp(:,i) = parts(ntotal_loc+nhalo_loc+i)%x(:)
        end do
        
        ! writing position
        call hdf5_parallel_write(fid,'virt/x',displ,global_dims,dbl_tmp)
        
        ! Converting velocity data to contiguous array
        do i = 1,nvirt_loc
            dbl_tmp(:,i) = parts(ntotal_loc+nhalo_loc+i)%vx(:)
        end do
        
        ! Writing velocity
        call hdf5_parallel_write(fid,'virt/v',displ,global_dims,dbl_tmp)
        
        ! cleanup
        deallocate(dbl_tmp)
        
        ! Writing data for ghost particles -----------------------------------------------------------------------------------------
        ! defining array shapes and displacments
        global_dims(:) = [dim,sum(nghos_glob)]
        displ(:) = [0,SUM(nghos_glob(1:procid))] ! hdf5 works with 0-indexing
        
        allocate(int_tmp(nghos_loc,1),dbl_tmp(nghos_loc,1))
        
        ! writing 1D arrays
        pos0 = ntotal_loc+nhalo_loc+nvirt_loc+1
        pos1 = ntotal_loc+nhalo_loc+nvirt_loc+nghos_loc
        int_tmp(:,1) = parts(pos0:pos1)%indglob
        call hdf5_parallel_write(fid,'ghos/ind',displ(2),global_dims(2),int_tmp(:,1))
        int_tmp(:,1) = procid
        call hdf5_parallel_write(fid,'ghos/procid',displ(2),global_dims(2),int_tmp(:,1))
        int_tmp(:,1) = parts(pos0:pos1)%itype
        call hdf5_parallel_write(fid,'ghos/type',displ(2),global_dims(2),int_tmp(:,1))
        dbl_tmp(:,1) = parts(pos0:pos1)%rho
        call hdf5_parallel_write(fid,'ghos/rho',displ(2),global_dims(2),dbl_tmp(:,1))
        dbl_tmp(:,1) = parts(pos0:pos1)%p
        call hdf5_parallel_write(fid,'ghos/p',displ(2),global_dims(2),dbl_tmp(:,1))
        
        deallocate( int_tmp,dbl_tmp )
        
        ! Converting position data to contiguous array
        allocate(dbl_tmp(dim,nghos_loc))
        do i = 1,nghos_loc
            dbl_tmp(:,i) = parts(ntotal_loc+nhalo_loc+nvirt_loc+i)%x(:)
        end do
        
        ! writing position
        call hdf5_parallel_write(fid,'ghos/x',displ,global_dims,dbl_tmp)
        
        ! Converting velocity data to contiguous array
        do i = 1,nghos_loc
            dbl_tmp(:,i) = parts(ntotal_loc+nhalo_loc+nvirt_loc+i)%vx(:)
        end do
        
        ! Writing velocity
        call hdf5_parallel_write(fid,'ghos/v',displ,global_dims,dbl_tmp)
        
        ! cleanup
        deallocate(dbl_tmp)
        
        allocate(dbl_tmp(tenselem,nghos_loc))
        
        global_dims(:) = [tenselem,sum(nghos_glob)]
        
        do i = 1,nghos_loc
            dbl_tmp(:,i) = parts(ntotal_loc+nhalo_loc+nvirt_loc+i)%sig(:)
        end do
        call hdf5_parallel_write(fid,'ghos/stress',displ,global_dims,dbl_tmp)
        
        deallocate(dbl_tmp)
        
        ! Closing output file
        call h5fclose_f(fid,ierr)
        
        ! closing hdf5
        call h5close_f(ierr)
        
	end subroutine output
	
	!==============================================================================================================================
	subroutine write_ini_config
	! Subroutine for writing initial configuration data using MPI IO subroutines
		
		implicit none
		integer:: ntotal_glob(numprocs)
		
		! Exchanging how many particles each process will output data for
		call MPI_ALLGATHER(ntotal_loc,1,MPI_INTEGER,ntotal_glob,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
		
        call h5open_f(ierr)
        
        filepath = trim(output_directory)//"/sph_out0000.h5"
        call hdf5_parallel_fileopen(filepath,fid)
        
        call h5gcreate_f(fid,"real",real_group_id,ierr)
        call h5gclose_f(real_group_id,ierr)
        
        call h5gcreate_f(fid,"halo",ghos_group_id,ierr)
        call h5gclose_f(ghos_group_id,ierr)
        
        call h5gcreate_f(fid,"virt",virt_group_id,ierr)
        call h5gclose_f(virt_group_id,ierr)
        
        call h5gcreate_f(fid,"ghos",ghos_group_id,ierr)
        call h5gclose_f(ghos_group_id,ierr)
        
        global_dims(:) = [dim,ntotal]

        allocate( int_tmp(ntotal_loc,1),dbl_tmp(ntotal_loc,1) )
        
        ! writing 1D arrays
        pos0 = 1
        pos1 = ntotal_loc
        int_tmp(:,1) = parts(pos0:pos1)%indglob
        call hdf5_parallel_write(fid,'real/ind',displ(2),global_dims(2),int_tmp(:,1))
        int_tmp(:,1) = procid
        call hdf5_parallel_write(fid,'real/procid',displ(2),global_dims(2),int_tmp(:,1))
        int_tmp(:,1) = parts(pos0:pos1)%itype
        call hdf5_parallel_write(fid,'real/type',displ(2),global_dims(2),int_tmp(:,1))
        dbl_tmp(:,1) = parts(pos0:pos1)%rho
        call hdf5_parallel_write(fid,'real/rho',displ(2),global_dims(2),dbl_tmp(:,1))
        dbl_tmp(:,1) = parts(pos0:pos1)%p
        call hdf5_parallel_write(fid,'real/p',displ(2),global_dims(2),dbl_tmp(:,1))
        
        deallocate( dbl_tmp,int_tmp )
        
        displ(:) = [0,SUM(ntotal_glob(1:procid))] ! hdf5 works with 0-indexing
        
        allocate(dbl_tmp(dim,ntotal_loc))
        do i = 1,ntotal_loc
            dbl_tmp(:,i) = parts(i)%x(:)
        end do
        call hdf5_parallel_write(fid,'real/x',displ,global_dims,dbl_tmp)
        do i = 1,ntotal_loc
            dbl_tmp(:,i) = parts(i)%vx(:)
        end do
        call hdf5_parallel_write(fid,'real/v',displ,global_dims,dbl_tmp)
        deallocate(dbl_tmp)
        
        ! Closing output file
        call h5fclose_f(fid,ierr)
        
        ! closing hdf5
        call h5close_f(ierr)
        
	end subroutine write_ini_config
	
end module output_m
