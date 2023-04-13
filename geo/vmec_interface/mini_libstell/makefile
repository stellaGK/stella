# Typically this makefile is only called from make in the parent directory,
# in which the variables FC and EXTRA_COMPILE_FLAGS are set.

LIBSTELL_OBJ_FILES = ezcdf.o ezcdf_inqvar.o handle_err.o stel_constants.o \
	ezcdf_GenGet.o ezcdf_opncls.o read_wout_mod.o safe_open_mod.o stel_kinds.o vsvd0.o \
	parse_extension.o

LIBSTELL_TARGET = mini_libstell.a

.PHONY: all clean

all: $(LIBSTELL_TARGET)

%.o: %.f90
	$(FC) $(EXTRA_COMPILE_FLAGS)  -DNETCDF -c $<

%.o: %.f
	$(FC) $(EXTRA_COMPILE_FLAGS) -DNETCDF -c $<

%.o: %.F
	$(FC) $(EXTRA_COMPILE_FLAGS) -DNETCDF -c $<

$(LIBSTELL_TARGET): $(LIBSTELL_OBJ_FILES)
	ar rcs $(LIBSTELL_TARGET) $(LIBSTELL_OBJ_FILES)

clean:
	rm -f *.o *.mod *.MOD *~ $(LIBSTELL_TARGET)

# Dependencies:
ezcdf.o: ezcdf_GenGet.o ezcdf_opncls.o
ezcdf_inqvar.o:
handle_err.o:
stel_constants.o: stel_kinds.o
ezcdf_GenGet.o: ezcdf_opncls.o ezcdf_inqvar.o
ezcdf_opncls.o:
read_wout_mod.o: stel_kinds.o safe_open_mod.o stel_constants.o vsvd0.o ezcdf.o parse_extension.o
safe_open_mod.o:
stel_kinds.o:
vsvd0.o: stel_kinds.o
parse_extension.o:
