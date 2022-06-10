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

FCFLAGS_gnu="-O3 -flto -fbacktrace -fimplicit-none -pedantic -Wall -Wextra -Jobj"
FCFLAGS_gnu_debug="-Og -g -fcheck=all -fbacktrace -fimplicit-none -pedantic -Wall -Wextra -Jobj"
FCFLAGS_intel="-O3 -ipo -traceback -module obj"
FCFLAGS_intel_debug="-g -traceback -checka all -module obj"
FCFLAGS_CUDA="-Mcuda -Minfo -Iobj -module obj"
FCFLAGS_CUDA_debug="-g -C -traceback -Iobj -module obj"
MKFILES_DIR=makefiles
MKFILE_serial=makefile.serial
MKFILE_MPI=makefile.MPI
MKFILE_CUDA=makefile.CUDA

case "$1" in
	"serial-gnu")
		FCFLAGS="$FCLAGS $FCFLAGS_gnu"
		MKFILE=$MKFILE_SERIAL
		;;
	"serial-intel")
		FCFLAGS="$FCFLAGS $FCFLAGS_intel"
		MKFILE=$MKFILE_SERIAL
		;;
	"serial-gnu-debug")
		FCFLAGS="$FCFLAGS $FCFLAGS_gnu_debug"
		MKFILE=$MKFILE_SERIAL
		;;
	"serial-intel-debug")
		FCFLAGS="$FCFLAGS $FCFLAGS_intel_debug"
		MKFILE=$MKFILE_SERIAL
		;;
	"MPI-gnu")
		FCFLAGS="$FCFLAGS $FCFLAGS_gnu"
		MKFILE=$MKFILE_MPI
		;;
	"MPI_intel")
		FCFLAGS="$FCFLAGS $FCFLAGS_intel"
		MKFILE=$MKFILE_MPI
		;;
	"MPI_gnu_debug")
		FCFLAGS="$FCFLAGS $FCFLAGS_gnu_debug"
		MKFILE=$MKFILE_MPI
		;;
	"MPI_intel_debug")
		FCFLAGS="$FCFLAGS $FCFLAGS_intel_debug"
		MKFILE=$mKFILE_MPI
		;;
	"CUDA")
		FCFLAGS="$FCFLAGS $FCFLAGS_CUDA"
		MKFILE=$MKFILE_CUDA
		;;
	"CUDA-debug")
		FCFLAGS="$FCFLAGS $FCFLAGS_CUDA_debug"
		MKFILE=$MKFILE_CUDA
		;;
	"clean")
		make -f $MKFILES_DIR/$MKFILE_serial clean
		exit 0
		;;
	*)
		echo "USAGE: $0 MODE"
	        echo "    where MODE = serial-gnu"
        	echo "                 serial-intel"
	        echo "                 serial-gnu-debug"
        	echo "                 serial-intel-debug"
	        echo "                 MPI-gnu"
        	echo "                 MPI-intel"
	        echo "                 MPI-gnu-debug"
        	echo "                 MPI-intel-debug"
        	echo "                 CUDA"
	        echo "                 CUDA-debug"
	        echo "relevant environment variables:"
        	echo "    FC: compiler to use"
	        echo "    FCFLAGS: flags for the compiler/linker"
		exit 1

esac

make -f $MKFILES_DIR/$MKFILE FCFLAGS="$FCFLAGS"