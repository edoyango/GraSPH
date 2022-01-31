FC=gfortran
SDIR=./src
ODIR=./obj
OPTIONS=-O3 -J$(ODIR)
_OBJ=constants.o \
param.o \
datatypes.o \
globvar.o \
summary_m.o \
material_rates_m.o \
output_m.o \
input_m.o \
kernel_m.o \
flink_lisk_m.o \
single_step_m.o \
time_integration_m.o \
main.o \

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.f90
	$(FC) -c -o $@ $< $(OPTIONS)

sph: $(OBJ)
	$(FC) -o $@ $^ $(OPTIONS)

clean:
	rm ./obj/*.o
	rm ./obj/*.obj
	rm ./obj/*.mod