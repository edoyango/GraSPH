#!/bin/bash
# script to compile sph binary.
# no arguments compile serial sph. One of mpi, gpu, omp compile

sdir_serial=./src_serial
sdir_mpi=./src_MPI
sdir_omp=./src_omp
sdir_gpu=./src_GPU
odir=./obj
outdir=./outputdata

# creating directories if needed
[ ! -d $odir ] && mkdir $odir
[ ! -d $outdir ] && mkdir $outdir

# checking how many "mode" arguments provided
mode=serial
modeCount=0
debug=0
clean=0
compiler=gnu
compilerCount=0
# loop to count number of compiler and mode args, and setting mode/compiler vars
for i in $@
do
    if [ $i == mpi ]; then
        let modeCount++
        mode=mpi
    fi
    if [ $i == gpu ]; then
        let modeCount++
        mode=gpu
    fi
    if [ $i == serial ]; then
        let modeCount++
        mode=serial
    fi
    if [ $i == gnu ]; then
        let compilerCount++
        compiler=gnu
    fi
    if [ $i == intel ]; then
        let compilerCount++
        compiler=intel
    fi
    if [ $i == nv ]; then
        let compilerCount++
        compiler=nv
    fi
    [ $i == clean ] && clean=1
    [ $i == debug ] && debug=1
done

# throwing error if conflicting args supplied
if [ $modeCount -gt 1 ]; then
    echo "Too many mode arguments specified! Specify only mpi, omp, gpu, or nothing for serial"
    exit 1
fi
if [ $compilerCount -gt 1 ]; then
    echo "Too many compiler arguments specified! Specify only intel, nv or nothing for gnu"
    exit 1
fi

# cleaning obj dir
if [ $clean == 1 ]; then
    echo "Clean argument supplied. Cleaning obj directory"
    rm $odir/*
    exit 0
fi

# setting compiler-appropriate compiler options and debug options if debug argument supplied
echo "Compilation mode: $mode"
if [ $mode == serial ] || [ $mode == mpi ]; then # serial/mpi mode (same flags)
    if [ $compiler == gnu ]; then
        [ $mode == serial ] && echo "compiler: gfortran"
        [ $mode == mpi ] && echo "compiler: mpifort (gfortran)"
        compoption="-O3 -flto -J$odir -I$odir"
        [ $debug == 1 ] && compoption="-Og -g -fcheck=all -fbacktrace -J$odir -I$odir"
    fi
    if [ $compiler == intel ]; then
        [ $mode == serial ] && echo "compiler: ifort"
        [ $mode == mpi ] && echo "compiler: mpifort (ifort)"
        compoption="-O3 -ipo -module $odir -I$odir"
        [ $debug == 1 ] && compoption="-g -traceback -check all -module $odir -I$odir"
    fi
    if [ $compiler == nv ]; then
        [ $mode == serial ] && echo "compiler: nvfortran"
        [ $mode == mpi ] && echo "compiler: mpifort (nvfortran)"
        compoption="-O3 -module $odir -I$odir"
        [ $debug ==1 ] && compoption="-g -C -traceback -module $odir -I$odir"
    fi
fi
if [ $mode == gpu ]; then # gpu mode
    echo "compiler: nvfortran"
    compoption="-module $odir -I$odir -Mcuda -Minfo"
    [ $debug == 1 ] && compoption="-g -C -traceback -module $odir -I$odir -Mcuda -Minfo"
fi
if [ $mode == omp ]; then # OpenMP mode
    if [ $compiler == gnu ]; then
        echo "compiler: gfortran"
        compoption="-O3 -flto -J$odir -I$odir -fopenmp"
        [ $debug == 1 ] && compoption="-Og -g -fcheck=all -fbacktrace -J$odir -I$odir -fopenmp"
    fi
    if [ $compiler == intel]; then
        echo "compiler: ifort"
        compoption="-O3 -ipo -module $odir -I$odir -openmp"
        [ $debug == 1 ] && compoption="-g -traceback -check all -module $odir -I$odir -openmp"
    fi
fi

echo "Flags: $compoption"
echo "Compile command:"
# Performing compilation based on mode and compiler
if [ $mode == serial ]; then
    srcfiles="$sdir_serial/param.f90 \
$sdir_serial/datatypes.f90 \
$sdir_serial/globvar.f90 \
$sdir_serial/summary_m.f90 \
$sdir_serial/output_m.f90 \
$sdir_serial/input_m.f90 \
$sdir_serial/kernel_m.f90 \
$sdir_serial/flink_lisk_m.f90 \
$sdir_serial/material_rates_m.f90 \
$sdir_serial/single_step_m.f90 \
$sdir_serial/time_integration_m.f90 \
$sdir_serial/main.f90"
    [ $compiler == gnu ] && FC=gfortran
    [ $compiler == intel ] && FC=ifort
    [ $compiler == nv ] && FC=nvfortran
fi
if [ $mode == mpi ]; then
    srcfiles="$sdir_mpi/param.f90 \
$sdir_mpi/datatypes.f90 \
$sdir_mpi/globvar.f90 \
$sdir_mpi/globvar_para.f90 \
$sdir_mpi/param_para.f90 \
$sdir_mpi/summary_m.f90 \
$sdir_mpi/error_msg_m.f90 \
$sdir_mpi/material_rates_m.f90 \
$sdir_mpi/output_m.f90 \
$sdir_mpi/input_m.f90 \
$sdir_mpi/kernel_m.f90 \
$sdir_mpi/flink_lisk_m.f90 \
$sdir_mpi/ORB_sr_m.f90 \
$sdir_mpi/ORB_m.f90 \
$sdir_mpi/single_step_m.f90 \
$sdir_mpi/time_integration_m.f90 \
$sdir_mpi/main.f90"
    FC=mpifort
fi
if [ $mode == gpu ]; then
    srcfiles="$sdir_gpu/param.cuf \
$sdir_gpu/datatypes.cuf \
$sdir_gpu/globvar.cuf \
$sdir_gpu/summary_m.cuf \
$sdir_gpu/output_m.cuf \
$sdir_gpu/input_m.cuf \
$sdir_gpu/kernel_m.cuf \
$sdir_gpu/flink_lisk_m.cuf \
$sdir_gpu/time_integration_m.cuf \
$sdir_gpu/main.cuf"
    FC=nvfortran
fi
echo "$FC $compoption $srcfiles -o sph"
$FC $compoption $srcfiles -o sph
