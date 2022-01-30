FC=mpifort
OPTIONS=-O3
SDIR=./src
ODIR=./obj
_OBJ=constants.o \
param.o \
datatypes.o \
globvar.o \
summary_m.o \
error_msg_m.o \
material_rates_m.o \
output_m.o \
input_m.o \
kernel_m.o \
flink_lisk_m.o \
ORB_sr_m.o \
ORB_m.o \
single_step_m.o \
time_integration_m.o \
main_sph.o \

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.f90
	$(FC) -c -o $@ $< $(OPTIONS)

sph: $(OBJ)
	$(FC) -o $@ $^ $(OPTIONS)

clean:
	rm ./obj/*.o
	rm ./obj/*.obj
