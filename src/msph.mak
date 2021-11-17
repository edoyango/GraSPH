FC=gfortran
OPTIONS= -O3
srcdir=
idir=./execs
bakdir=./execs.bak
objects=main_sph.o \
	art_visc.o \
	density.o \
	dsearch.o \
	ext_force.o \
	input.o \
	int_force.o \
	integration.o \
	kernel.o \
	output.o \
	single_step.o \
	time_print.o \
	virt_part.f90
#
mpsh: $(objects)
	$(FC) $(OPTIONS) -o mcgsph.exe $(objects)
#
	if [ -d $(bakdir) ]; then \
	echo "execs.bak Directory Exists"; else \
	mkdir $(bakdir); \
	echo "execs.bak Directory Created"; \
	fi
#
	if [ -d $(idir) ]; then \
	echo "execs Directory Exists"; else \
	mkdir $(idir); \
	echo "execs Directory Created"; \
	fi
#
clean:
	rm *.o
	rm *.obj
	rm *~
#
%.o: %.f90
	$(FC) $(OPTIONS) -c -o $@ $<
