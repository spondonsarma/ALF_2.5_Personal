.PHONY : all

all: cov_eq.out cov_scal.out cov_tau.out  cov_tau_ph.out \
     Max_SAC.out Max_SAC_T0.out \

cov_tau_ph.f90 : cov_tau.f90
	cpp -w -P -DPartHole cov_tau.f90 -o cov_tau_ph.f90

%.out: %.f90
	$(ALF_FC) -c $(ALF_FLAGS_ANA) $<
	$(ALF_FC) -o $@ $*.o $(ALF_LIB)

clean:
	rm -f *.mod *~ \#* *.o   cov_tau_ph.f90

clean_all:
	rm -f *.mod *~ \#* *.o *.out  cov_tau_ph.f90
