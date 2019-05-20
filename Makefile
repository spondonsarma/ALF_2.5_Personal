.PHONY : all lib ana program
#Hub_Ising SPT Hub Hub_Can Kondo_Honey
all: lib ana program
#Hub_Ising SPT  Hub_Can Kondo_Honey Z2_Slave

lib:
	cd Libraries && $(MAKE)
ana: lib
	cd Analysis && $(MAKE)
program: lib
	cd Prog && $(MAKE) Examples
#Hub_Ising: lib
#	cd Prog && $(MAKE) Hub_Ising
#SPT: lib
#	cd Prog && $(MAKE) SPT
#Hub_Can: lib
#	cd Prog && $(MAKE) Hub_Can
#Kondo_Honey: lib
#	cd Prog && $(MAKE) Kondo_Honey
Z2_Slave: lib
	cd Prog && $(MAKE) Z2_Slave

.PHONY : clean cleanall cleanprog cleanlib cleanana help
clean: cleanall
cleanall: cleanprog cleanlib cleanana  
cleanprog:
	cd Prog && $(MAKE) clean 
cleanlib:
	cd Libraries && $(MAKE) clean
cleanana:
	cd Analysis && $(MAKE) clean
help:
	@echo "The following are some of the valid targets of this Makefile"
	@echo "all, program, lib, ana, clean, cleanall, cleanprog, cleanlib, cleanana"
#	@echo "Hub_Ising SPT Hub Hub_Can Kondo_Honey"

