OBJS=zgeqr2.o zgeqrf.o zlarf.o zlarfb.o zlarfg.o zlarft.o zung2r.o zungqr.o

all: $(OBJS)
	ar -r libqrref.a $(OBJS)

.SUFFIXES: .f90
.f90.o:
	$(FC) $(SUFFIX) $(FLAGS) $< 

clean: 
	rm $(OBJS)
