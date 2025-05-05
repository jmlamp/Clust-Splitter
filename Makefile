
#==============================
# Makefile for Clust-Splitter
#==============================

# Requires a Fortran compiler (here gfortran) to be installed.
#
# After running this makefile, run the program by typing: ajo.exe

FF = gfortran -O2

FFLAGS = -fbounds-check -Wall -fcheck=all -g -ffpe-trap=invalid,zero,overflow -fbacktrace -pg

OPEN = 

RM = del   #in windows
#RM = rm    #in unix/linux/mac

all: ajo

ajo: clustsplitter.o
	$(FF) -o ajo $(FFLAGS) $(OPEN) clustsplitter.o
	
clustsplitter.o: clustsplitter.f95
	$(FF) -c $(FFLAGS) $(OPEN) clustsplitter.f95