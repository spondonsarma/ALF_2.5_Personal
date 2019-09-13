.PHONY : Compile tidy clean Examples Z2_Slave Z2_Matter Hub_Can Kondo_Honey

OBJS= LRC_mod.o control_mod.o Fields_mod.o Operator.o WaveFunction.o observables_mod.o Predefined_structures.o \
      $(HAMILTONIAN) QDRP_decompose.o udv_state.o Hop_mod.o UDV_WRAP.o Global_mod.o Wrapgr_mod.o \
      tau_m.o tau_p.o main.o wrapul.o cgr1.o wrapur.o cgr2_2.o upgrade.o 
MODS= control.mod fields_mod.mod global_mod.mod hop_mod.mod lrc_mod.mod observables.mod \
      operator_mod.mod predefined_structures.mod qdrp_mod.mod tau_m_mod.mod tau_p_mod.mod \
      udv_state_mod.mod udv_wrap_mod.mod wavefunction_mod.mod wrapgr_mod.mod hamiltonian.mod
OHAM= Hamiltonian_Examples.o Hamiltonians/Hamiltonian_Z2_slave_spins.o Hamiltonians/Hamiltonian_Z2_Matter.o \
      Hamiltonians/Hamiltonian_Hub_Canonical.o Hamiltonians/Hamiltonian_Kondo_Honey.o
BINS= Examples.out Z2_slave_spins.out Z2_Matter.out Hubb_Can.out Kondo_Honey.out
      
all:
	(make tidy Examples tidy Z2_Slave tidy Z2_Matter tidy Hub_Can tidy Kondo_Honey)

Examples:
	./Git_config.sh
	(make Compile TARGET="Examples.out" HAMILTONIAN="Hamiltonian_Examples.o")

Z2_Slave:
	./Git_config.sh
	(make Compile TARGET="Z2_slave_spins.out" HAMILTONIAN="Hamiltonians/Hamiltonian_Z2_slave_spins.o")

Z2_Matter:
	./Git_config.sh
	(make Compile TARGET="Z2_Matter.out" HAMILTONIAN="Hamiltonians/Hamiltonian_Z2_Matter.o")

Hub_Can:
	./Git_config.sh
	(make Compile TARGET="Hubb_Can.out" HAMILTONIAN="Hamiltonians/Hamiltonian_Hub_Canonical.o")

Kondo_Honey:
	./Git_config.sh
	(make Compile TARGET="Kondo_Honey.out" HAMILTONIAN="Hamiltonians/Hamiltonian_Kondo_Honey.o")


Compile: $(OBJS)
	$(ALF_FC) -o $(TARGET) $(OBJS) $(ALF_LIB)

.SUFFIXES: .F90 .f
.f.o .F90.o:
	$(ALF_FC) -c -o $@ $(ALF_FLAGS_PROG) $<

tidy:
	rm -f $(OBJS) $(MODS) $(OHAM)

clean: tidy
	rm -f $(BINS)
