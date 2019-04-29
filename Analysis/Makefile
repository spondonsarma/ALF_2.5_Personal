.PHONY : all
all: cov_eq.out cov_scal.out cov_tau.out jackv5.out \
     Max_SAC.out Max_SAC_ph.out Max_SAC_T0.out \

%.out: %.f90
	$(ALF_FC) -c $(ALF_FLAGS_ANA) $<
	$(ALF_FC) -o $@ $*.o $(ALF_LIB)

clean:
	rm -f *.mod *~ \#* *.o
