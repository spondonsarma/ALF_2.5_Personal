.PHONY : all Examples Z2_Slave Hub_Ising Hub_Can SPT Kondo_Honey NematicDirac clean
all: Examples Z2_Slave Hub_Ising Hub_Can SPT Kondo_Honey NematicDirac

Examples:
	./Git_config.sh
	make -f Compile_Examples

Z2_Slave:
	./Git_config.sh
	make -f Compile_Z2_Slave

Hub_Ising:
	./Git_config.sh
	make -f Compile_Hub_Ising

Hub_Can:
	./Git_config.sh
	make -f Compile_Hub_Can

SPT:
	./Git_config.sh
	make -f Compile_SPT

Kondo_Honey:
	./Git_config.sh
	make -f Compile_Kondo_Honey

NematicDirac:
	./Git_config.sh
	make -f Compile_NematicDirac


clean: 	
	(make -f Compile_Examples   clean )
	(make -f Compile_Hub_Can   clean )
	(make -f Compile_Hub_Ising   clean )
	(make -f Compile_SPT   clean )
	(make -f Compile_Kondo_Honey clean )
	(make -f Compile_Z2_Slave clean )
	(make -f Compile_NematicDirac clean )
	rm -f *.mod *~ \#* 
