.PHONY : all tidy clean
BINS=cov_eq.out cov_scal.out cov_tau.out cov_tau_ph.out Max_SAC.out ana.out #repack_latt.out

all: $(OBJS1) $(BINS)

cov_tau_ph.F90 : cov_tau.F90
	cpp -w -P -DPartHole cov_tau.F90 -o cov_tau_ph.F90

%.o: %.F90
	$(ALF_FC) -c -o $@ $(ALF_FLAGS_ANA) $<

%.out: Predefined_Latt_mod.o ana_mod.o %.o
	$(ALF_FC) -o $@ $^ $(ALF_LIB)

tidy:
	rm -f ana_mod.mod predefined_lattices.mod cov_tau_ph.F90

clean: tidy
	rm -f $(BINS)
