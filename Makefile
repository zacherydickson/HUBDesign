CC = g++

ALN = mafft
LCR = dustmasker

all: SA_BOND/src/sabond.cpp config
	${CC} -fopenmp SA_BOND/src/sabond.cpp -o SA_BOND/sa_bond

config: bin/Config.pl
	bin/Config.pl -a ${ALN} -l ${LCR} > HUBDesign.cfg

clean:
	@echo "Cleaning up binaries and configuration files"
	@rm -rvf SA_BOND/sa_bond *.cfg
