
HOST=$(shell hostname)

ifeq ($(HOST),RSFAS-22676.local)
MATLAB_HOME=/Applications/MATLAB_R2016b.app/
MATLAB_LIB = $(MATLAB_HOME)/bin/maci64/
endif

ifeq ($(HOST),RSFAS-22673)
MATLAB_HOME=/home/honglin/local/MATLAB/R2016b
MATLAB_LIB = $(MATLAB_HOME)/bin/glnxa64/

endif

# To compile and run, you need to run the following command in shell
# Desktop: export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/honglin/local/MATLAB/R2016b/bin/glnxa64/
# MAC    : #export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Applications/MATLAB_R2016b.app/bin/maci64/
a.out: main1.cpp BOUNDLP.hpp
	g++ -std=c++11 -g -I$(MATLAB_HOME)/extern/include/ -L$(MATLAB_LIB) $^ -lmx -lmat 
