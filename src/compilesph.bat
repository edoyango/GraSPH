ifort /traceback /O3 /Qxhost /Qipo param.f90 arrayAlloc.f90 art_visc.f90 density.f90 ext_force.f90 flink_lisk.f90 input.f90 int_force.f90 integration.f90 kernel.f90 main_sph.f90 output.f90 single_step.f90 time_print.f90 virt_part.f90 /exe:runsph.exe

del *.obj *.mod