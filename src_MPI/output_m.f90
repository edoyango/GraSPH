module output_m

	use globvar,		only: ntotal_loc,nhalo_loc,nvirt_loc,parts,itimestep,save_step
	use globvar_para,	only: procid,numprocs,ierr,MPI_ftype
	use mpi
	use param,			only: f,dim,output_directory,output_phys,output_halo,output_virt,output_flt_type
	
	implicit none
	character(len=220),private:: filepath
	character(len=4),private::number
	integer,allocatable,private:: int_tmp(:,:)
	real(f),allocatable,private:: dbl_tmp(:,:)
	integer(kind=MPI_OFFSET_KIND),private:: displ_int,displ_dbl
	integer,private:: n,i,j,k,d,fh,request(3),pos0,pos1,status(3,MPI_STATUS_SIZE),intsize,dblsize
	
	public:: output,write_ini_config
	private:: write_data
	
contains

	!==============================================================================================================================
	subroutine output
	! Subroutine to write data to disk. Called intermittently, determined by save_step supplied at run-time
		
		implicit none
		integer:: ntotal_glob(numprocs),nhalo_glob(numprocs),nvirt_glob(numprocs)
		
		! Exchanging how many particles each process will output data for
		call MPI_IALLGATHER(ntotal_loc,1,MPI_INTEGER,ntotal_glob,1,MPI_INTEGER,MPI_COMM_WORLD,request(1),ierr)
		call MPI_IALLGATHER(nhalo_loc,1,MPI_INTEGER,nhalo_glob,1,MPI_INTEGER,MPI_COMM_WORLD,request(2),ierr)
		call MPI_IALLGATHER(nvirt_loc,1,MPI_INTEGER,nvirt_glob,1,MPI_INTEGER,MPI_COMM_WORLD,request(3),ierr)
		
		!Format number
		n = itimestep/save_step
		write(number,1000) n
		1000 format(I4.4)
		
		! intializing integer, double sizes
		intsize = SIZEOF(1); dblsize = SIZEOF(1_f)
		
		call MPI_WAITALL(3,request,status,ierr)
		
		displ_int = SUM(ntotal_glob(1:procid))*intsize
		displ_dbl = SUM(ntotal_glob(1:procid))*dblsize
		
		! Writing misc particle info: ID, process owned by, itype -----------------------------------------------------------------
		allocate( int_tmp(3,ntotal_loc) )
		do i = 1,ntotal_loc
			int_tmp(1,i) = parts(i)%indglob
			int_tmp(2,i) = procid
			int_tmp(3,i) = parts(i)%itype
		end do
		
		filepath = trim(output_directory)//"/f_ind"//number//".dat"
		call write_data(int_tmp,dbl_tmp,displ_int,TRIM(filepath),LEN(TRIM(FILEPATH)),'int')
		
		! Writing particle position and velocity ----------------------------------------------------------------------------------
		if (output_phys(1)) then
			
			allocate( dbl_tmp(2*dim,ntotal_loc) )
			do i = 1,ntotal_loc
				dbl_tmp(1:dim,i) = parts(i)%x(1:dim)
				dbl_tmp(dim+1:2*dim,i) = parts(i)%vx(1:dim)
			end do
			
			filepath = trim(output_directory)//"/f_xv"//number//".dat"
			call write_data(int_tmp,dbl_tmp,displ_dbl,TRIM(filepath),LEN(TRIM(FILEPATH)),output_flt_type)
		
		end if
		
		! Writing density, pressure -----------------------------------------------------------------------------------------------
		if (output_phys(2)) then
			
			allocate( dbl_tmp(2,ntotal_loc) )
			do i = 1,ntotal_loc
				dbl_tmp(1,i) = parts(i)%rho
				dbl_tmp(2,i) = parts(i)%p
			end do
			
			filepath = trim(output_directory)//"/f_state"//number//".dat"
			call write_data(int_tmp,dbl_tmp,displ_dbl,TRIM(filepath),LEN(TRIM(FILEPATH)),output_flt_type)
		
		end if
		
		displ_int = SUM(nhalo_glob(1:procid))*intsize
		displ_dbl = SUM(nhalo_glob(1:procid))*dblsize
		
		if (output_halo(1)) then
			
			pos0 = ntotal_loc + 1
			pos1 = ntotal_loc + nhalo_loc
			
			! integer stuff: global particle index, procid, itype
			allocate( int_tmp(3,nhalo_loc) )
			do i = pos0,pos1
				int_tmp(1,i-ntotal_loc) = parts(i)%indglob
				int_tmp(2,i-ntotal_loc) = procid
				int_tmp(3,i-ntotal_loc) = parts(i)%itype
			end do
			
			filepath = trim(output_directory)//"/h_ind"//number//".dat"
			call write_data(int_tmp,dbl_tmp,displ_int,TRIM(filepath),LEN(TRIM(FILEPATH)),'int')
			
			! Position, velocity
			allocate( dbl_tmp(2*dim,nhalo_loc) )
			do i = pos0,pos1
				dbl_tmp(1:dim,i-ntotal_loc) = parts(i)%x(1:dim)
				dbl_tmp(dim+1:2*dim,i-ntotal_loc) = parts(i)%vx(1:dim)
			end do
			
			filepath = trim(output_directory)//"/h_xv"//number//".dat"
			call write_data(int_tmp,dbl_tmp,displ_dbl,TRIM(filepath),LEN(TRIM(FILEPATH)),output_flt_type)
			
			! density, pressure
			if (output_halo(2)) then
				
				allocate( dbl_tmp(2,nhalo_loc) )
				do i = pos0,pos1
					dbl_tmp(1,i-ntotal_loc) = parts(i)%rho
					dbl_tmp(2,i-ntotal_loc) = parts(i)%p
				end do
				
				filepath = trim(output_directory)//"/h_state"//number//".dat"
				call write_data(int_tmp,dbl_tmp,displ_dbl,TRIM(filepath),LEN(TRIM(FILEPATH)),output_flt_type)
				
			end if
			
		end if
		
		! Writing virtual particle data -------------------------------------------------------------------------------------------
		displ_int = SUM(nvirt_glob(1:procid))*intsize
		displ_dbl = SUM(nvirt_glob(1:procid))*dblsize
		if (output_virt(1)) then
			
			pos0 = ntotal_loc + nhalo_loc + 1
			pos1 = ntotal_loc + nhalo_loc + nvirt_loc
			
			! integer stuff: global particle index, procid, itype
			allocate( int_tmp(3,nvirt_loc) )
			do i = pos0,pos1
				int_tmp(1,i-ntotal_loc-nhalo_loc) = parts(i)%indglob
				int_tmp(2,i-ntotal_loc-nhalo_loc) = procid
				int_tmp(3,i-ntotal_loc-nhalo_loc) = parts(i)%itype
			end do
			
			filepath = trim(output_directory)//"/v_ind"//number//".dat"
			call write_data(int_tmp,dbl_tmp,displ_int,TRIM(filepath),LEN(TRIM(FILEPATH)),'int')
			
			! Position, velocity
			allocate( dbl_tmp(2*dim,nvirt_loc) )
			do i = pos0,pos1
				dbl_tmp(1:dim,i-ntotal_loc-nhalo_loc) = parts(i)%x(1:dim)
				dbl_tmp(dim+1:2*dim,i-ntotal_loc-nhalo_loc) = parts(i)%vx(1:dim)
			end do
			
			filepath = trim(trim(output_directory)//"/v_xv"//number//".dat")
			call write_data(int_tmp,dbl_tmp,displ_dbl,TRIM(filepath),LEN(TRIM(FILEPATH)),output_flt_type)
			
			! density, pressure
			if (output_virt(2)) then
				
				allocate( dbl_tmp(2,nvirt_loc) )
				do i = pos0,pos1
					dbl_tmp(1,i-ntotal_loc-nhalo_loc) = parts(i)%rho
					dbl_tmp(2,i-ntotal_loc-nhalo_loc) = parts(i)%p
				end do
				
				filepath = trim(output_directory)//"/v_state"//number//".dat"
				call write_data(int_tmp,dbl_tmp,displ_dbl,TRIM(filepath),LEN(TRIM(FILEPATH)),output_flt_type)
				
			end if
			
		end if
		
	end subroutine output
	
	!==============================================================================================================================
	subroutine write_ini_config
	! Subroutine for writing initial configuration data using MPI IO subroutines
		
		implicit none
		integer:: ntotal_glob(numprocs)
		
		! Exchanging how many particles each process will output data for
		call MPI_ALLGATHER(ntotal_loc,1,MPI_INTEGER,ntotal_glob,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
		
		displ_int = SUM(ntotal_glob(1:procid))*SIZEOF(1)
		displ_dbl = SUM(ntotal_glob(1:procid))*SIZEOF(1_f)
		
		! particle index, process, itype data -------------------------------------------------------------------------------------
		allocate( int_tmp(3,ntotal_loc) )
		do i = 1,ntotal_loc
			int_tmp(1,i) = parts(i)%indglob
			int_tmp(2,i) = procid
			int_tmp(3,i) = parts(i)%itype
		end do
		
		filepath = trim(output_directory)//"/ini_ind.dat"
		call write_data(int_tmp,dbl_tmp,displ_int,TRIM(filepath),LEN(TRIM(FILEPATH)),'int')
		
		! particle position and velocity data -------------------------------------------------------------------------------------
		allocate( dbl_tmp(2*dim,ntotal_loc) )
		do i = 1,ntotal_loc
			dbl_tmp(1:dim,i) = parts(i)%x(1:dim)
			dbl_tmp(dim+1:2*dim,i) = parts(i)%vx(1:dim)
		end do
		
		filepath = trim(output_directory)//"/ini_xv.dat"
		call write_data(int_tmp,dbl_tmp,displ_dbl,TRIM(filepath),LEN(TRIM(FILEPATH)),output_flt_type)
		
		! density and pressure ----------------------------------------------------------------------------------------------------
		allocate( dbl_tmp(2,ntotal_loc) )
		
		do i = 1,ntotal_loc
			dbl_tmp(1,i) = parts(i)%rho
			dbl_tmp(2,i) = parts(i)%p
		end do
		
		filepath = trim(output_directory)//"/ini_state.dat"
		call write_data(int_tmp,dbl_tmp,displ_dbl,TRIM(filepath),LEN(TRIM(FILEPATH)),output_flt_type)
		
	end subroutine write_ini_config
	
	! Subroutine to write dbl data arrays -----------------------------------------------------------------------------------------
	subroutine write_data(int_data,flt_data,npart_displ,filepath,nchar,datatype)
	! Helper subroutine to write data. Can write as double or single precision floats, or integer data
		
		implicit none
		integer,intent(in):: nchar
		integer,allocatable,intent(inout):: int_data(:,:)
		real(f),allocatable,intent(inout):: flt_data(:,:)
		integer(kind=MPI_OFFSET_KIND),intent(in):: npart_displ
		character(len=3),intent(in):: datatype
		integer(kind=MPI_OFFSET_KIND):: displ
		integer:: dims(2),fh
		character(len=nchar):: filepath
		
		call MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(filepath),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
		
		select case (datatype)
			case('dbl')
				dims = SHAPE(flt_data(:,:))
				displ = dims(1)*npart_displ
				call MPI_FILE_WRITE_AT_ALL(fh,displ,flt_data,dims(1)*dims(2),MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
			case('sgl')
				dims = SHAPE(flt_data(:,:))
				displ = dims(1)*npart_displ
				call MPI_FILE_WRITE_AT_ALL(fh,displ,flt_data,dims(1)*dims(2),MPI_REAL,MPI_STATUS_IGNORE,ierr)
			case('int')
				dims = SHAPE(int_data(:,:))
				displ = dims(1)*npart_displ
				call MPI_FILE_WRITE_AT_ALL(fh,displ,int_data,dims(1)*dims(2),MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
		end select
		
		call MPI_FILE_CLOSE(fh,ierr)
		
		if (allocated(int_data)) deallocate(int_data)
		if (allocated(flt_data)) deallocate(flt_data)
	
	end subroutine write_data
	
end module output_m
