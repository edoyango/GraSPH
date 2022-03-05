@echo off

@echo beginning compile script

REM setting up variables
set serial_src_dir=.\src_serial\
set MPI_src_dir=.\src_MPI\
set OMP_src_dir=.\src_OMP\
set obj_dir=.\obj\
set out_dir=.\outputdata\
set compoption=/O3 /Qipo /traceback

REM creating directories if needed
if not exist %obj_dir% (
	@echo Obj directory does not exist. Creating...
	mkdir %obj_dir% )
if not exist %out_dir% (
	@echo Outputdata directory does not exist. Creating...
	mkdir %out_dir% )

REM checking how many "mode" arguments provided
set mode=serial
set modeCount=0
set debug=0
for %%x in (%*) do (
	if "%%x"=="mpi" (
		set mode=mpi
		set /A modeCount+=1 )
	if "%%x"=="omp" (
		set mode=omp
		set /A modeCount+=1 )
	if "%%x"=="gpu" (
		set mode=gpu
		set /A modeCount+=1 )
	if "%%x"=="clean" (
		@echo Cleaning object directory
		goto clean )
	if "%%x"=="debug" (
		set debug=1 )
)

REM throw error if too many modes supplied
if %modeCount% gtr 1 (
	@echo Too many arguments specified! Specify only mpi, omp, gpu, or nothing for serial
	exit /b )
	
REM set compile options based on mode
if %modeCount%==1 (
	REM throw error if gpu mode is supplied
	if %mode%==gpu (
		@echo gpu argument supplied. Cannot compile gpu code - HPC SDK and therefore nvfortran only available on Linux
		exit /b )
	@echo Compiling SPH code in %mode% mode
	if %mode%==omp set compoption=%compoption% /Qopenmp 
	if %mode%==mpi set compoption=%compoption% /Qxhost 
)

REM compiling for serial code if no mode supplied
if %modeCount%==0 (
	@echo Compiling SPH code in %mode% mode
	set compoption=%compoption% /Qxhost )
if %debug%==1 (
	set compoption=/check:all /traceback )

@echo Compiling with %compoption%

if %mode%==serial goto:serial
if %mode%==mpi goto:mpi
if %mode%==omp goto:omp

REM section for compilation of serial code
:serial
	@echo Compiling source code in src_serial...
	ifort /object:%obj_dir% /module:%obj_dir% %compoption%^
	%serial_src_dir%constants.f90 ^
	%serial_src_dir%param.f90 ^
	%serial_src_dir%datatypes.f90 ^
	%serial_src_dir%globvar.f90 ^
	%serial_src_dir%summary_m.f90 ^
	%serial_src_dir%output_m.f90 ^
	%serial_src_dir%input_m.f90 ^
	%serial_src_dir%kernel_m.f90 ^
	%serial_src_dir%flink_lisk_m.f90 ^
	%serial_src_dir%material_rates_m.f90 ^
	%serial_src_dir%single_step_m.f90 ^
	%serial_src_dir%time_integration_m.f90 ^
	%serial_src_dir%main.f90 /exe:.\sph.exe
	goto:eof

REM section for compilation of MPI code
:mpi
	@echo Compiling source code in src_MPI...
	mpiifort /object:%obj_dir% /module:%obj_dir% %compoption%^
	%MPI_src_dir%param.f90 ^
	%MPI_src_dir%datatypes.f90 ^
	%MPI_src_dir%globvar.f90 ^
	%MPI_src_dir%globvar_para.f90 ^
	%MPI_src_dir%param_para.f90 ^
	%MPI_src_dir%summary_m.f90 ^
	%MPI_src_dir%error_msg_m.f90 ^
	%MPI_src_dir%material_rates_m.f90 ^
	%MPI_src_dir%output_m.f90 ^
	%MPI_src_dir%input_m.f90 ^
	%MPI_src_dir%kernel_m.f90 ^
	%MPI_src_dir%flink_lisk_m.f90 ^
	%MPI_src_dir%ORB_sr_m.f90 ^
	%MPI_src_dir%ORB_m.f90 ^
	%MPI_src_dir%single_step_m.f90 ^
	%MPI_src_dir%time_integration_m.f90 ^
	%MPI_src_dir%main.f90 /exe:.\sph.exe
	goto:eof

REM section for compilation of OMP code
:omp
	@echo Compiling source code in src_omp...
	ifort /object:%obj_dir% /module:%obj_dir% %compoption%^
	%OMP_src_dir%constants.f90 ^
	%OMP_src_dir%param.f90 ^
	%OMP_src_dir%datatypes.f90 ^
	%OMP_src_dir%globvar.f90 ^
	%OMP_src_dir%summary_m.f90 ^
	%OMP_src_dir%output_m.f90 ^
	%OMP_src_dir%input_m.f90 ^
	%OMP_src_dir%kernel_m.f90 ^
	%OMP_src_dir%flink_lisk_m.f90 ^
	%OMP_src_dir%material_rates_m.f90 ^
	%OMP_src_dir%single_step_m.f90 ^
	%OMP_src_dir%time_integration_m.f90 ^
	%OMP_src_dir%main.f90 /exe:.\sph.exe
	goto:eof 

:clean
del %obj_dir%\*.obj
del %obj_dir%\*.mod
del %obj_dir%\*.o