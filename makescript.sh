#!/bin/bash
# this script is a layer on top of a collection of make files.
# It handles the logic of choosing the variant to compile.
# first argument is variant which is one of:
#	serial-gnu
#	serial-intel
#	serial-gnu-debug
#	serial-intel-debug
#	MPI-gnu
#	MPI-intel
#	MPI-gnu-debug
#	MPI-intel-debug
#	CUDA
#	CUDA-debug

FCFLAGS_gnu_opt="-O3 -flto"
FCFLAGS_gnu_debug="-Og -g -fcheck=all"
FCFLAGS_gnu="-fbacktrace -fimplicit-none -pedantic -Wall -Wextra -Jobj"
FCFLAGS_intel_opt="-O3 -ipo"
FCFLAGS_intel_debug="-g -check all"
FCFLAGS_intel="-traceback -module obj -warn all"
FCFLAGS_CUDA_opt=""
FCFLAGS_CUDA_debug="-g -C"
FCFLAGS_CUDA="-Mcuda -traceback -Iobj -module obj -Minfo"
MKFILES_DIR=makefiles
MKFILE_serial=makefile.serial
MKFILE_MPI=makefile.mpi
MKFILE_CUDA=makefile.cuda

make -f $MKFILES_DIR/$MKFILE_serial clean

SERIAL=1
MPI=0
CUDA=0
GNU=1
OPTIMISE=0
DEBUG=0
INTEL=0
for arg in $@
do
	[ $arg = "--mpi" ] && MPI=1 && SERIAL=0
	[ $arg = "--cuda" ] && CUDA=1 && SERIAL=0 && GNU=0
	[ $arg = "--optimise" ] && OPTIMISE=1
	[ $arg = "--debug" ] && DEBUG=1
	[ $arg = "--intel" ] && GNU=0 && INTEL=1
	if [ $arg = "--help" ] || [ $arg = "-h" ]; then
		echo "Running this script without any options will compile the Fortran code in src_serial into"
		echo " the executable sph-serial. The following options will augment the behavour:"
		echo "--mpi      : compile the parallel code in src_MPI into the executable sph-mpi."
		echo "--cuda     : compile the CUDA Fortran coda in src_GPU into the executable sph-cuda."
		echo "--optimise : compile with optimisation options. Which options depend on which compiler and"
		echo "             source code being compiled."
		echo "--debug    : compile with debug options. Which options depend on which compiler and source"
		echo "             code being compiled."
	        echo "relevant environment variables:"
        	echo "    FC: compiler to use"
	        echo "    FCFLAGS: flags for the compiler/linker"
		exit 1
	fi
done

[ $MPI -eq 1 ] && [ $CUDA -eq 1 ] && \
	echo "Only --cuda or --mpi (or neither) should be supplied!" && return 1
[ $INTEL -eq 1 ] && [ $CUDA -eq 1 ] && \
	echo "Only --cuda or --intel (or neither) should be supplied!" && return 1

[ $SERIAL -eq 1 ] && MKFILE=$MKFILE_serial
[ $MPI -eq 1 ] && MKFILE=$MKFILE_MPI
[ $CUDA -eq 1 ] && MKFILE=$MKFILE_CUDA

if [ $INTEL -eq 1 ]; then
	FCFLAGS="$FCFLAGS $FCFLAGS_intel"
	[ $OPTIMISE -eq 1 ] && FCFLAGS="$FCFLAGS $FCFLAGS_intel_opt"
	[ $DEBUG -eq 1 ] && FCFLAGS="$FCFLAGS $FCFLAGS_intel_debug"
elif [ $GNU -eq 1 ]; then
	FCFLAGS="$FCFLAGS $FCFLAGS_gnu"
	[ $OPTIMISE -eq 1 ] && FCFLAGS="$FCFLAGS $FCFLAGS_gnu_opt"
	[ $DEBUG -eq 1 ] && FCFLAGS="$FCFLAGS $FCFLAGS_gnu_debug"
elif [ $CUDA -eq 1 ]; then
	FCFLAGS="$FCFLAGS $FCFLAGS_CUDA"
	echo $FCFLAGS
	[ $OPTIMISE -eq 1 ] && FCFLAGS="$FCFLAGS $FCFLAGS_CUDA_opt"
	[ $DEBUG -eq 1 ] && FCFLAGS="$FCFLAGS $FCFLAGS_CUDA_debug"
fi


echo $MKFILES_DIR/$MKFILE
if [ -v FC ]
then
	make -f $MKFILES_DIR/$MKFILE FCFLAGS="$FCFLAGS" FC="$FC"
else
	make -f $MKFILES_DIR/$MKFILE FCFLAGS="$FCFLAGS"
fi
