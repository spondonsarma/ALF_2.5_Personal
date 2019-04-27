LIB=modules_90.a
OBJS=mat_mod.o Random_Wrap.o errors.o Files_mod.o maxent.o  matrix.o  maxent_stoch.o  fourier.o \
         Histogram.o  lattices_v3.o  Natural_constants.o log_mesh.o  precdef.mod.o \
         Histogram_v2.o Mat_subroutines.o

$(LIB): $(OBJS)
	ar  -r $(LIB) $(OBJS)

%.o: %.f90
	$(ALF_FC) -c $(ALF_FLAGS_MODULES) $<

%.o: %.F90
	$(ALF_FC) -c $(ALF_FLAGS_MODULES) $<

errors.o: mat_mod.o

maxent.o: mat_mod.o errors.o

lattices_v3.o: matrix.o

clean:
	rm -f *.mod $(OBJS) *~ \#* 

