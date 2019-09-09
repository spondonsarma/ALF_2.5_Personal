.PHONY : all Examples Z2_Slave Hub_Ising Hub_Can SPT Kondo_Honey NematicDirac clean
all: Examples Z2_Slave Hub_Ising Hub_Can SPT Kondo_Honey NematicDirac

OBJS= LRC_mod.o control_mod.o Fields_mod.o Operator.o WaveFunction.o observables_mod.o Predefined_structures.o \
      $(HAMILTONIAN) QDRP_decompose.o udv_state.o Hop_mod.o UDV_WRAP.o Global_mod.o Wrapgr_mod.o \
      tau_m.o tau_p.o main.o wrapul.o cgr1.o wrapur.o cgr2_2.o upgrade.o 


Examples:
	./Git_config.sh
	(make Compile TARGET="Examples.out" HAMILTONIAN="Hamiltonian_Examples.o")
# 	make -f Compile_Examples

Z2_Slave:
	./Git_config.sh
	(make Compile TARGET="Hamiltonian_Z2_slave_spins.out" HAMILTONIAN="Hamiltonians/Hamiltonian_Z2_slave_spins.o")
# 	make -f Compile_Z2_Slave

Z2_Matter:
	./Git_config.sh
	(make Compile TARGET="Hamiltonian_Z2_Matter.out" HAMILTONIAN="Hamiltonians/Hamiltonian_Z2_Matter.o")

Hub_Can:
	./Git_config.sh
	(make Compile TARGET="Hubb_Can.out" HAMILTONIAN="Hamiltonians/Hamiltonian_Hub_Canonical.o")
# 	make -f Compile_Hub_Can

SPT:
	./Git_config.sh
	(make Compile TARGET="SPT.out" HAMILTONIAN="Hamiltonians/Hamiltonian_SPT.o")
# 	make -f Compile_SPT

Kondo_Honey:
	./Git_config.sh
	(make Compile TARGET="Kondo_Honey.out" HAMILTONIAN="Hamiltonians/Hamiltonian_Kondo_Honey_mod.o")
# 	make -f Compile_Kondo_Honey

NematicDirac:
	./Git_config.sh
	(make Compile TARGET="NematicDirac.out" HAMILTONIAN="Hamiltonians/Hamiltonian_NematicDirac.o")
# 	make -f Compile_NematicDirac


Compile: $(OBJS)
	$(ALF_FC) -o $(TARGET) $(OBJS) $(ALF_LIB)

.SUFFIXES: .f90 .f .F90
.f.o .f90.o .F90.o:
	$(ALF_FC) -c -o $@ $(ALF_FLAGS_PROG) $<

clean:
	rm -f $(OBJS) *.mod *.smod *.o *~ \#* Hamiltonians/*.o

# clean: 	
# 	(make -f Compile_Examples   clean )
# 	(make -f Compile_Hub_Can   clean )
# 	(make -f Compile_Hub_Ising   clean )
# 	(make -f Compile_SPT   clean )
# 	(make -f Compile_Kondo_Honey clean )
# 	(make -f Compile_Z2_Slave clean )
# 	(make -f Compile_NematicDirac clean )
# 	rm -f *.mod *.o *~ \#* Hamiltonians/*.o
