.PHONY : Compile tidy clean Z2_Matter Kondo Hubbard tV Hubbard_Plain_Vanilla LRC .FORCE

OBJS= Hamiltonians/LRC_mod.o  Set_random_mod.o control_mod.o Fields_mod.o Operator_mod.o WaveFunction_mod.o observables_mod.o ContainerElementBase_mod.o DynamicMatrixArray_mod.o OpTTypes_mod.o\
      Predefined_Int_mod.o Predefined_Obs_mod.o Predefined_Latt_mod.o Predefined_Hop_mod.o  Predefined_Trial_mod.o \
      $(HAMILTONIAN) QDRP_decompose_mod.o udv_state_mod.o Hop_mod.o UDV_WRAP_mod.o \
      wrapul_mod.o cgr1_mod.o Global_mod.o wrapur_mod.o Langevin_HMC_mod.o upgrade_mod.o Wrapgr_mod.o cgr2_2_mod.o tau_m_mod.o tau_p_mod.o main.o 
MODS= control.mod fields_mod.mod global_mod.mod hop_mod.mod lrc_mod.mod observables.mod \
      operator_mod.mod predefined_int.mod predefined_lattices.mod predefined_obs.mod \
      predefined_hoppings.mod  predefined_trial.mod qdrp_mod.mod tau_m_mod.mod tau_p_mod.mod \
      udv_state_mod.mod udv_wrap_mod.mod wavefunction_mod.mod wrapgr_mod.mod hamiltonian.mod \
      containerelementbase_mod.mod opttypes_mod.mod dynamicmatrixarray_mod.mod \
			cgr1_mod.mod cgr2_2_mod.mod langevin_hmc_mod.mod set_random_mod.mod upgrade_mod.mod \
			wrapul_mod.mod wrapur_mod.mod
OHAM= Hamiltonians/Hamiltonian_Z2_Matter_mod.o \
	Hamiltonians/Hamiltonian_Kondo_mod.o \
	Hamiltonians/Hamiltonian_Hubbard_mod.o \
	Hamiltonians/Hamiltonian_Hubbard_Plain_Vanilla_mod.o \
	Hamiltonians/Hamiltonian_tV_mod.o \
	Hamiltonians/Hamiltonian_LRC_mod.o
# Add your new Hamiltonian here by appending  Hamiltonians/Hamiltonian_My_New_Hamiltonian_mod.o to the above list.
BINS= Z2_Matter.out Kondo.out Hubbard.out Hubbard_Plain_Vanilla.out tV.out \
      LRC.out
# Add your new Hamiltonian by appending My_New_Hamiltonian.out to the above list
      
all: $(BINS)

#My_New_Hamiltonian:
#	(make Compile TARGET="My_New_Hamiltonian.out" HAMILTONIAN="My_New_Hamiltonian_mod.o")
#Add you new Hamiltonian by uncommenting the above

Hubbard: Hubbard.out
Hubbard_Plain_Vanilla: Hubbard_Plain_Vanilla.out
tV: tV.out
Kondo: Kondo.out
LRC: LRC.out
Z2_Matter: Z2_Matter.out

%.out: .FORCE
	(make Compile TARGET="$*.out" HAMILTONIAN="Hamiltonians/Hamiltonian_$*_mod.o")

Compile: $(OBJS)
	$(ALF_FC) -o $(TARGET) $(OBJS) $(ALF_LIB)

%.o: %.F90
	$(ALF_FC) -c -o $@ $(ALF_FLAGS_PROG) $<

tidy:
	rm -f $(OBJS) $(MODS) $(OHAM)

clean: tidy
	rm -f $(BINS)

git.h: .FORCE
	./Git_config.sh

$(HAMILTONIAN): .FORCE Fields_mod.o Hamiltonians/LRC_mod.o observables_mod.o \
	Operator_mod.o Predefined_Hop_mod.o Predefined_Int_mod.o Predefined_Latt_mod.o \
	Predefined_Obs_mod.o Predefined_Trial_mod.o WaveFunction_mod.o

Hamiltonians/LRC_mod.o:
cgr1_mod.o: QDRP_decompose_mod.o udv_state_mod.o UDV_WRAP_mod.o
cgr2_2_mod.o: QDRP_decompose_mod.o udv_state_mod.o UDV_WRAP_mod.o
ContainerElementBase_mod.o: 
control_mod.o: 
DynamicMatrixArray_mod.o: ContainerElementBase_mod.o
Fields_mod.o: 
Global_mod.o: cgr1_mod.o control_mod.o Fields_mod.o observables_mod.o \
	Operator_mod.o udv_state_mod.o UDV_WRAP_mod.o wrapul_mod.o $(HAMILTONIAN)
Hop_mod.o: ContainerElementBase_mod.o DynamicMatrixArray_mod.o OpTTypes_mod.o \
	Operator_mod.o $(HAMILTONIAN)
Langevin_HMC_mod.o: cgr1_mod.o control_mod.o Hop_mod.o udv_state_mod.o \
	wrapul_mod.o wrapur_mod.o $(HAMILTONIAN)
main.o: git.h cgr1_mod.o control_mod.o Fields_mod.o Global_mod.o Hop_mod.o \
	Langevin_HMC_mod.o Operator_mod.o Set_random_mod.o tau_m_mod.o tau_p_mod.o \
	udv_state_mod.o Wrapgr_mod.o wrapul_mod.o wrapur_mod.o $(HAMILTONIAN)
observables_mod.o: 
Operator_mod.o: Fields_mod.o
OpTTypes_mod.o: ContainerElementBase_mod.o Operator_mod.o
Predefined_Hop_mod.o: Operator_mod.o WaveFunction_mod.o
Predefined_Int_mod.o: Operator_mod.o
Predefined_Latt_mod.o: 
Predefined_Obs_mod.o: observables_mod.o
Predefined_Trial_mod.o: Operator_mod.o Predefined_Hop_mod.o WaveFunction_mod.o
QDRP_decompose_mod.o: 
Set_random_mod.o: 
tau_m_mod.o: cgr2_2_mod.o control_mod.o Hop_mod.o Langevin_HMC_mod.o \
	Operator_mod.o udv_state_mod.o wrapur_mod.o $(HAMILTONIAN)
tau_p_mod.o: cgr1_mod.o control_mod.o Hop_mod.o Langevin_HMC_mod.o Operator_mod.o \
	tau_m_mod.o udv_state_mod.o wrapur_mod.o $(HAMILTONIAN)
udv_state_mod.o: QDRP_decompose_mod.o
UDV_WRAP_mod.o: QDRP_decompose_mod.o
upgrade_mod.o: control_mod.o Fields_mod.o $(HAMILTONIAN)
WaveFunction_mod.o: 
Wrapgr_mod.o: control_mod.o Fields_mod.o Hop_mod.o Operator_mod.o upgrade_mod.o \
	$(HAMILTONIAN)
wrapul_mod.o: Hop_mod.o udv_state_mod.o UDV_WRAP_mod.o $(HAMILTONIAN)
wrapur_mod.o: Hop_mod.o udv_state_mod.o $(HAMILTONIAN)
