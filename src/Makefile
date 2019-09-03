#OpenMM_INSTALL_DIR=/home/ymwang/apps/openmm_src/openmm_bin
#OpenMM_INSTALL_DIR=/export/apps/OpenMM/7.2.1
#OpenMM_INSTALL_DIR=/home/ymwang/apps/anaconda2/pkgs/openmm-7.1.1-py27_0
OpenMM_INSTALL_DIR=/home/ymwang/apps/anaconda2/pkgs/openmm-7.2.1-py27_1
ROOT=$(shell pwd)

CFLAGS = -g -std=c++11 -O3
#CXX=mpicxx
CXX=mpiicpc
LIB_DIR=$(OpenMM_INSTALL_DIR)/lib
INCLUDE_DIR=$(OpenMM_INSTALL_DIR)/include

LIBS= -lOpenMM
ALL_HEADERS = *.h
kernel.out: $(ALL_HEADERS) main.cpp
	$(CXX) $(CFLAGS) -I$(INCLUDE_DIR) main.cpp -L$(LIB_DIR) $(LIBS) -o kernel.out
