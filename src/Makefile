#
# Makefile for Quantum Quartic Hamiltonian Problem
#
# by Currix TM
#
BINPATH = ../bin/
####################################################################
quartic_HAM_HO_SRC = nrtype.f90 defparam_quartic_HO.f90 quartic_potential_HO.f90
quartic_HAM_HO_OBJ = $(quartic_HAM_HO_SRC:.f90=.o)
####################################################################
#
.SUFFIXES: .o .f90 .oOMP 
###
FC = gfortran
FOPT	= -O3 -Wall
PREP = -cpp
#
MODLIB   = -I/usr/local/lib/lapack95/lapack95_modules -I/usr/local/include/lapack95_modules
LAPACK95 = -L/usr/local/lib/lapack95 -llapack95 -L/usr/local/lib -llapack95
LAPACK77 = -llapack
BLAS     = -lblas
#
LIBS     = $(LAPACK95)  $(LAPACK77) $(BLAS)
# 
LIB_INTLIB = -lintlib
######################################################################################################
all: quartic_HO_energies
#################################################################
.PHONY : all
#################################################################
.f90.o:
	$(info )
	$(info Compiling single thread object file:)
	$(FC) -c $(PREP) $(FOPT) $(MODLIB) -o "$@" "$<"
#################################################################
.f90.oOMP:
	$(info )
	$(info Compiling multi thread object file:)
	$(FC) -c -fopenmp $(PREP) $(MODLIB) -o "$@" "$<"
#################################################################
#######################
quartic_HO_energies: $(quartic_HAM_HO_OBJ) Makefile 
	$(FC) $(FOPT) $(MODLIB) -o $(BINPATH)/$@_$(FC) $(quartic_HAM_HO_OBJ) $(LIBS) 
#######################
#######################
clean:
	$(info Cleaning object, module, and executable files.)
	@rm -f  *.o *.mod $(BINPATH)/quartic_HO_energies_$(FC)
#######################




