set src_dir=./src/
set obj_dir=./obj/
ifort ^
/object:%obj_dir% /module:%obj_dir% ^
/traceback /O3 /Qxhost /Qipo ^
%src_dir%constants.f90 ^
%src_dir%param.f90 ^
%src_dir%datatypes.f90 ^
%src_dir%globvar.f90 ^
%src_dir%arrayAlloc.f90 ^
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
%src_dir%time_print.f90 ^
%src_dir%virt_part.f90 ^
/exe:./runsph.exe