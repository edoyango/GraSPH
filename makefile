SHELL=/bin/bash -O extglob -c
FC=nvfortran# compiler
SDIR=./src_GPU# source file directory
ODIR=./obj# directory to store .o, .mod, .obj files generated by compilation

OPTIONS=-module $(ODIR) -I$(ODIR) -Mcuda -Minfo#-Og -J$(ODIR) -flto
_OBJ=param.o \
datatypes.o \
globvar.o \
summary_m.o \
output_m.o \
input_m.o \
kernel_m.o \
flink_lisk_m.o \
time_integration_m.o \
main.o \

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

sph: $(OBJ)
	$(FC) -o $@ $^ $(OPTIONS)

clean:
	@rm $(ODIR)/*(*.o|*.obj|*.mod)

$(ODIR)/%.o: $(SDIR)/%.cuf
	@if [ ! -d $(ODIR) ]; then \
		mkdir $(ODIR); \
	fi
	@if [ ! -d outputdata ]; then \
		mkdir outputdata; \
	fi
	$(FC) -c -o $@ $< $(OPTIONS)
