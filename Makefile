.PHONY : all lib ana program test_pgi_bug test_compiler_set

all: lib ana program

lib: test_pgi_bug
	cd Libraries && $(MAKE)
ana: lib
	cd Analysis && $(MAKE)
program: lib
	cd Prog && $(MAKE) program

.PHONY : clean cleanall cleanprog cleanlib cleanana tidy tidyprog tidyana help
clean: cleanall
cleanall: cleanprog cleanlib cleanana
tidy: cleanlib tidyana tidyprog
cleanprog:
	cd Prog && $(MAKE) clean 
cleanlib:
	cd Libraries && $(MAKE) clean
cleanana:
	cd Analysis && $(MAKE) clean
tidyana:
	cd Analysis && $(MAKE) tidy
tidyprog:
	cd Prog && $(MAKE) tidy

help:
	@echo "The following are some of the valid targets of this Makefile"
	@echo "all, lib, ana, program, clean, cleanall, cleanprog, cleanlib, cleanana"

test_compiler_set:
	@if [ -z ${ALF_FC} ]; then \
          printf "\n\033[0;31m Environment variable ALF_FC not set.\n"; \
          printf " Please source configure.sh before compilation.\033[0m\n"; \
          exit 1; \
        fi
	@if [ ! `command -v ${ALF_FC}` ]; then \
	  printf "\n\033[0;31m Compiler '${ALF_FC}' does not exist.\n"; \
	  printf " Please install or choose other in configure.sh.\033[0m\n"; \
	  exit 1; \
	fi

test_pgi_bug: test_compiler_set pgi_bug.F90
	@printf "Testing for specific bug in older PGI compilers..."
	@$(ALF_FC) -o pgi_bug.out pgi_bug.F90
	@./pgi_bug.out || \
	  (printf "\n\033[0;31mPGI bug present, aborting compilation.\033[0m\n"; \
	  exit 1; )
	@rm -f pgi_bug.out test_mod.mod
	@printf "done\n"
