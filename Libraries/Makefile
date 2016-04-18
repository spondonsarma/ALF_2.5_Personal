FC= $(mpif90)
FC= $(f90)
FLAGS= $(FL)
LIBS= $(Libs)/Modules/modules_90.a \
      $(Libs)/MyEis/libeis.a \
      $(Libs)/MyNag/libnag.a \
      $(Libs)/MyLin/liblin.a \
      $(LIB_BLAS_LAPACK)

all:
	(cd Modules;make FC=$(FC) FLAGS="$(FLAGS)") ;\
	(cd MyEis;make FC=$(FC) FLAGS="$(FLAGS)") ;\
	(cd MyNag;make FC=$(FC) FLAGS="$(FLAGS)") ;\
	(cd MyLin;make FC=$(FC) FLAGS="$(FLAGS)")

clean:
	(cd Modules;make clean);\
	(cd MyEis;make clean);\
	(cd MyNag;make clean);\
	(cd MyLin;make clean)

