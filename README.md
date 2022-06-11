# GraSPH
GraSPH is an SPH program originally inteded for simulations of bulk granular material as well as fluids. 
This repo contains Fortran source code files in src_serial, src_MPI, src_GPU, seperated by the level of paralellism.
- src_serial contains code intended only on one core. Compilation requires only compiler (confirmed working with gfortran 9.4.0, ifort 2021.6.0).
- src_MPI contains code intended to run multi-core configuration enabled with the MPI communication standard. (confirmed working with gfortran 9.4.0, ifort 2021.6.0, openmpi 4.1.2)
- src_GPU contains code intended to run on a CUDA-enabled GPU.

# Prerequisites
Compilation of all source code requires hdf5 libraries where the serial and cuda code requires the high-level libraries as well.
- hdf5 (with high-level libraries and built with openmpi)
- make
- An appropriate compiler:
  - nvfortran (tested with v22.1.0) - mandatory for GPU code;
  - gfortran (tested with v9.4.0); or
  - ifort (tested with v2021.6.0).
- An MPI implementation (for the MPI code only)

# Compiling
Compiling the code is done via the makescript.sh script. This can be invoked with one of the following arguments:
- `makescript.sh serial-gnu`: compiles the code  in src_serial into the executable `sph-serial` using optimisation options specific to gfortran
- `makescript.sh serial-intel`: compiles the code in src_serial into the executable `sph-serial` using optimisation options specific to ifort
- `makescript.sh mpi-gnu`: compiles the code in src_MPI into the executable `sph-mpi` using optimisation options specific to gfortran (requires an mpi implementation)
- `makescript.sh mpi-intel`: compiles the code in src_MPI into the executable `sph-mpi` using optimisation options specific to ifort (requires an mpi implementation)
- `makescript.sh cuda`: compiles the code in src_GPU into the executable `sph-cuda` using optimisation options specific to nvfortran (requires NVIDIA HPC SDK).

# Running
Running the program requires running the executable and supplying three integer arguments e.g.:
`./sph-serial <max timesteps> <print interval> <write interval>`
where <max timesteps> is the maximum number of timesteps to run the simulation for, 
  <print interval> is the number of steps between printing to the terminal, and
  <write interval> is the number of steps between writing data to disk.
The first compilation should run the classic dam-break experiment with 65,000 SPH particles. This can be run after first compilation by
`./sph-serial 7500 100 100`

 Currently geometry and simulation parameters are hardcoded. These can be controlled by
 - virt_part subroutine in input.f90 (boundary geometry)
 - input subroutine in input.f90 (initial geometry)
 - param.f90 (parameters).
 Work is ongoing to improve this.
 
# Visualisation
Currently visualisation is done via a MatLab script which parses the output hdf5 files. This is included as Plot_hdf5.m. This workflow shall be improved in the future.
