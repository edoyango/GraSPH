module output_m

    use globvar,        only: ntotal_loc,nhalo_loc,nvirt_loc,nghos_loc,parts,itimestep,save_step,ntotal,nvirt,maxnloc
    use globvar_para,    only: procid,numprocs,ierr,MPI_ftype
    use param,            only: f,dim,output_directory,output_phys,output_halo,output_virt,output_flt_type
    
    use hdf5
    use mpi
    
    use hdf5_parallel_io_helper_m, only: hdf5_parallel_fileopen,hdf5_parallel_write
    
    implicit none
    character(len=220),private:: filepath
    character(len=4),private::number
    integer,private:: n,i,request(4),status(MPI_STATUS_SIZE),posrange(2)
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
        
        call h5gcreate_f(fid,"virt",virt_group_id,ierr)
        call h5gclose_f(virt_group_id,ierr)
        
        call h5gcreate_f(fid,"halo",ghos_group_id,ierr)
        call h5gclose_f(ghos_group_id,ierr)
        
        call h5gcreate_f(fid,"ghos",ghos_group_id,ierr)
        call h5gclose_f(ghos_group_id,ierr)
        
        ! Writing data for real particles ------------------------------------------------------------------------------------------
        ! defining array shapes and displacments
        call MPI_WAIT(request(1),status,ierr) ! Waiting for non-blocking exchange to complete
        
        global_dims(:) = [dim,ntotal]
        displ(:) = [0,SUM(ntotal_glob(1:procid))] ! hdf5 works with 0-indexing
        
        posrange = [1,ntotal_loc]
        call write_particle_data(fid,'real',posrange,global_dims,displ)
        
        ! Writing data for halo particles ------------------------------------------------------------------------------------------
        ! defining array shapes and displacments
        call MPI_WAIT(request(2),status,ierr) ! Waiting for non-blocking exchange to complete
        
        global_dims(:) = [dim,sum(nhalo_glob)]
        displ(:) = [0,SUM(nhalo_glob(1:procid))] ! hdf5 works with 0-indexing
        
        posrange(:) = ntotal_loc + [1,nhalo_loc]
        call write_particle_data(fid,'halo',posrange,global_dims,displ)
        
        ! Writing data for virtual particles ---------------------------------------------------------------------------------------
        ! defining array shapes and displacments
        call MPI_WAIT(request(3),status,ierr) ! Waiting for non-blocking exchange to complete
        
        global_dims(:) = [dim,sum(nvirt_glob)]
        displ(:) = [0,SUM(nvirt_glob(1:procid))] ! hdf5 works with 0-indexing
        
        posrange(:) = ntotal_loc + nhalo_loc + [1,nvirt_loc]
        call write_particle_data(fid,'virt',posrange,global_dims,displ)

        ! Writing data for ghost particles -----------------------------------------------------------------------------------------
        ! defining array shapes and displacments
        call MPI_WAIT(request(4),status,ierr) ! Waiting for non-blocking exchange to complete
        
        global_dims(:) = [dim,sum(nghos_glob)]
        displ(:) = [0,SUM(nghos_glob(1:procid))] ! hdf5 works with 0-indexing
        
        posrange(:) = ntotal_loc + nhalo_loc + nvirt_loc + [1,nghos_loc]
        call write_particle_data(fid,'ghos',posrange,global_dims,displ)
        
        ! Closing output file
        call h5fclose_f(fid,ierr)
        
        ! closing hdf5
        call h5close_f(ierr)
        
    end subroutine output
        
    !==============================================================================================================================
    subroutine write_ini_config
    ! Subroutine for writing initial configuration data using MPI IO subroutines
        
        implicit none
        integer:: ntotal_glob(numprocs),posrange(2)
        
        
        ! Exchanging how many particles each process will output data for
        call MPI_IALLGATHER(ntotal_loc,1,MPI_INTEGER,ntotal_glob,1,MPI_INTEGER,MPI_COMM_WORLD,request(1),ierr)
        
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
        
        call MPI_WAIT(request(1),status,ierr)
        
        global_dims(:) = [dim,ntotal]
        displ(:) = [0,SUM(ntotal_glob(1:procid))] ! hdf5 works with 0-indexing
        
        posrange(:) = [1,ntotal_loc]
        call write_particle_data(fid,'real',posrange,global_dims,displ)
        
        ! Closing output file
        call h5fclose_f(fid,ierr)
        
        ! closing hdf5
        call h5close_f(ierr)
        
    end subroutine write_ini_config

    !==============================================================================================================================
    subroutine write_particle_data(fid_in,particletype,posrange_in,gdims,ldispl)
    
        implicit none
        integer(HID_T),intent(in):: fid_in
        character(*),intent(in):: particletype
        integer,intent(in):: posrange_in(2)
        integer(HSIZE_T),intent(in):: gdims(2)
        integer(HSSIZE_T),intent(in):: ldispl(2)
        integer:: nelem
        integer,allocatable:: int_tmp(:,:)
        real(f),allocatable:: dbl_tmp(:,:)
        
        nelem = posrange_in(2)-posrange_in(1)+1
        allocate(int_tmp(nelem,1),dbl_tmp(nelem,1))
        
        ! writing 1D arrays
        int_tmp(:,1) = parts(posrange_in(1):posrange_in(2))%indglob
        call hdf5_parallel_write(fid_in,particletype//'/ind',ldispl(2),gdims(2),int_tmp(:,1))
        int_tmp(:,1) = procid
        call hdf5_parallel_write(fid_in,particletype//'/procid',ldispl(2),gdims(2),int_tmp(:,1))
        int_tmp(:,1) = parts(posrange_in(1):posrange_in(2))%itype
        call hdf5_parallel_write(fid_in,particletype//'/type',ldispl(2),gdims(2),int_tmp(:,1))
        dbl_tmp(:,1) = parts(posrange_in(1):posrange_in(2))%rho
        call hdf5_parallel_write(fid_in,particletype//'/rho',ldispl(2),gdims(2),dbl_tmp(:,1))
        dbl_tmp(:,1) = parts(posrange_in(1):posrange_in(2))%p
        call hdf5_parallel_write(fid_in,particletype//'/p',ldispl(2),gdims(2),dbl_tmp(:,1))
        
        deallocate(int_tmp,dbl_tmp)
        
        ! Converting position data to contiguous array
        allocate(dbl_tmp(dim,nelem))
        
        do i = 1,nelem
            dbl_tmp(:,i) = parts(posrange_in(1)+i-1)%x(:)
        end do
        
        ! writing position
        call hdf5_parallel_write(fid,particletype//'/x',ldispl,gdims,dbl_tmp)
        
        ! Converting velocity data to contiguous array
        do i = 1,nelem
            dbl_tmp(:,i) = parts(posrange_in(1)+i-1)%vx(:)
        end do
        
        ! writing velocity
        call hdf5_parallel_write(fid,particletype//'/v',ldispl,gdims,dbl_tmp)
        
        ! cleanup
        deallocate(dbl_tmp)
        
    end subroutine write_particle_data

end module output_m
