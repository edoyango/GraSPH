# GraSPH
GraSPH is an SPH program originally inteded for simulations of bulk granular material as well as fluids. 
This repo contains Fortran source code files in src_CAF and src_GPU, seperated by the level of paralellism. This code is an upgraded version of used in \[[1](https://doi.org/10.1016/j.compgeo.2020.103474), [2](https://doi.org/10.1007/s11440-021-01162-4), [3](https://doi.org/10.26180/14484099.v1)\]. Which incorporates more of Fortran's features, namely derived types and coarrays. It also incorporates structural changes enabling faster run times, and an option for GPU acceleration. The code currently only simulates water via the classic "Weakly-compressible" approach. Simple granular models will be implemented soon.
- src_CAF contains code intended to run multi-core configuration enabled with the Coarray Fortran 2008 features. (confirmed working with gfortran 9.4.0, ifort 2021.8.0). Serial runs can also be setup (explained later).
- src_GPU contains code intended to run on a CUDA-enabled GPU.

NOTE This is a hobby project and is actively being developed. Major changes to the main branch can occur.

# Prerequisites
Compilation of all source code requires hdf5 libraries where the serial and cuda code requires the high-level libraries as well.
- hdf5 (with high-level libraries and built with openmpi)
- make
- An appropriate compiler:
  - nvfortran (tested with v21.9.0 to 23.0.0) - mandatory for GPU code;
  - gfortran (tested with v9.4.0 to v11.2.0); or
  - ifort (tested with v2021.8.0)/ifx (tested with v2023.0.0).
- An MPI implementation (for the MPI code only)
- If using `gfortran`, [OpenCoarrays](https://github.com/sourceryinstitute/OpenCoarrays) library is needed.

# Compiling
Compiling the code is done via the makefiles in the `makefiles` directory. This can be invoked with one of the following arguments:

## CoArray Fortran
the compilation command is
```
make FC=<compiler> compiler=<gnu/intel> mode=<serial/caf> extras=<opt,dev,debug,singlecaf> -f makefiles/makefile.caf
```
where
`FC` specifies the compiler to use e.g., gfortran, ifort, caf, mpifort, or mpiifort etc. (default is gfortran)
`compiler` specifies which type of compiler i.e., either gnu or intel (default is gnu).
`mode` specifies whether you would like to compile in serial or using Coarrays. 
`extras` add extra options:
* `opt` adds further optimisation options (which reduce usefulness of error messages)
* `dev` adds compiler warnings
* `debug` adds more detailed error messages and checks (which reduces performance)
* `singlecaf` is useful only when compiling with `mode=intel`. It compiles the executable with `-coarray=single`, but uses parallel IO i.e., needs to link to HDF5 parallel IO libraries. This allows the `sph` executable to be executed with `mpiexec`.

## CUDA Fortran
the compilation command is
```
make -f makefiles/makefile.cuda
```
This already assumes `nvfortran` as the compiler. Currently doesn't allow extra customisation (unless you modify the makefile). This will soon be merged with the CAF makefile.

## Important Environment Variables
`FCFLAGS` is used to specify compiler options. You should set this to point to Fortran modulefiles. e.g., `export FCFLAGS="-I/usr/include/hdf5/openmpi"` points the compiler to the HDF5 OpenMPI Fortran modules installed with `apt-get`.
`LDFLAGS` is used to specify linker options. You should set this to point to shared libraries. e.g., `export LDFLAGS="-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi"` points the linker to the HDF5 OpenMPI shared libraries installed by `apt-get` and needed by GraSPH. Note that `-lhdf5` is included by default, but you may also need to add options specifying the HL libraries (lhdf5_hl_fortran.so or lhdf5hl_fortran.so), in the case of compiling with `mode=serial`) or the HDF5_fortran libraries (libhdf5_fortran.so).

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
 
## CAF (`mode=caf`)
Built with Opencoarrays
`cafrun -n <ncpus> ./sph 100000 1000 1000`

Built with `ifort`/`ifx`:
`FOR_COARRAY_NUM_IMAGES=<ncpus> ./sph 100000 1000 1000`
or if built with `extras=singlecaf`
`mpiexec -n <ncpus> ./sph 100000 1000 1000`

# Example Build and Run (with gfortran and Opencoarrays)
Installing opencoarrays and openmpi via spack (brew-linux is also a good option)
```
$ spack install opencoarrays
```
Making sure `caf` and `cafrun` are in my path
```
$ spack load opencoarrays openmpi
$ which caf
/usr/local/spack/opt/spack/linux-ubuntu20.04-skylake/gcc-9.4.0/opencoarrays-2.7.1-wiecvev57rcwa6wdobdhmk2fukurcm6d/bin/caf
$ which cafrun
/usr/local/spack/opt/spack/linux-ubuntu20.04-skylake/gcc-9.4.0/opencoarrays-2.7.1-wiecvev57rcwa6wdobdhmk2fukurcm6d/bin/cafrun
```
If the above `which` commands don't return anything, you'll need to add the caf bin directory into your path.

Compiling the program
```
$ cd ~/GraSPH
# Compiling using caf (FC=caf), with coarrays active (mode=caf), with extra optimisation and development warnings information (extras=opt,dev)
$ make FC=caf mode=caf extras=opt,dev -f makefiles/makefile.caf
```
Running the program
```
$ cafrun -n 4 sph 100000 1000 1000
```

# Visualisation
Currently visualisation is done via a MatLab script which parses the output hdf5 files. This is included as Plot_hdf5.m. This workflow shall be improved in the future.
