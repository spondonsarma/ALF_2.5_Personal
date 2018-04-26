.PHONY : all
all: jackv5.out cov_tau.out cov_scal.out Max_SAC.out \
     Max_SAC_ph.out cov_eq.out 

%.out: %.f90
	$(ALF_FC) -o $*.o $(ALF_FLAGS_ANA) $<
	$(ALF_FC) -o $@ $*.o $(ALF_LIB)

clean:
	rm -f *.mod *~ \#* *.out *.o
