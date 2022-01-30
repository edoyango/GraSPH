@echo off

set src_dir=.\src\
set obj_dir=.\obj\

REM checking if debug is passed to script if so, compile with /check option
set compoption=/traceback /O3 /Qxhost /Qipo
:loop
	if "%1"=="debug" (
		set compoption=/traceback /O3 /Qxhost /Qipo /check:all )
	if "%1"=="clean" (
		goto clean )
	if "%1"=="defaultoption" (
		set compoption=/traceback )
shift
if not "%~1"=="" goto loop

ifort ^
/object:%obj_dir% /module:%obj_dir% %compoption% ^
%src_dir%constants.f90 ^
%src_dir%param.f90 ^
%src_dir%datatypes.f90 ^
%src_dir%globvar.f90 ^
%src_dir%summary_m.f90 ^
%src_dir%art_visc.f90 ^
%src_dir%density.f90 ^
%src_dir%ext_force.f90 ^
%src_dir%flink_lisk.f90 ^
%src_dir%input.f90 ^
%src_dir%int_force.f90 ^
%src_dir%integration.f90 ^
%src_dir%kernel.f90 ^
%src_dir%main_sph.f90 ^
%src_dir%output.f90 ^
%src_dir%single_step.f90 ^
%src_dir%virt_part.f90 ^
/exe:./runsph.exe

goto:eof

:clean
del %obj_dir%\*.obj
del %obj_dir%\*.mod