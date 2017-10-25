#****************************
#*** Copyright Notice ***
#IMPACT-Z,Copyright (c) 2016, The Regents of the University of California, through
#Lawrence Berkeley National Laboratory (subject to receipt of any required approvals
#from the U.S. Dept. of Energy).  All rights reserved.
#If you have questions about your rights to use or distribute this software,
#please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.
#NOTICE.  This Software was developed under funding from the U.S. Department of Energy
#and the U.S. Government consequently retains certain rights. As such, the U.S. Government
#has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable,
#worldwide license in the Software to reproduce, distribute copies to the public, prepare
#derivative works, and perform publicly and display publicly, and to permit other to do so.
#****************************

#*****************************************************
#  General Makefile
#
#*****************************************************

#**************************************************************************
# Macros defining the Fortran, C/C++ compiler and linker.

CC = mpif90
LINK = mpif90
FFLAGS = -O3

#**************************************************************************
# List of .o files that EXENAME depends on.  Edit as appropriate for MP.

OBJS = \
        TPSA/mod_mathfunc.o TPSA/mod_polymap.o TPSA/mod_tpsa.o \
	DataStruct/NumConst.o DataStruct/PhysConst.o DataStruct/Pgrid.o \
	DataStruct/Data.o \
	Func/Timer.o  \
	Appl/DriftTube.o Appl/ConstFoc.o \
	Appl/BeamLineElem.o \
	Appl/CompDom.o  Appl/BeamBunch.o  Appl/Distribution.o \
	Contrl/Input.o Contrl/Output.o Contrl/AccSimulator.o Contrl/main.o

OBJS2 = \
	mod_mathfunc.o mod_polymap.o mod_tpsa.o \
	NumConst.o PhysConst.o Pgrid.o Data.o \
        Timer.o \
	DriftTube.o ConstFoc.o \
	BeamLineElem.o CompDom.o BeamBunch.o Distribution.o \
	Input.o Output.o AccSimulator.o main.o	


#**************************************************************************
# Change this line if you don't like 'a.out'.

EXEOLD = ImpactZexe

#************************************************************************
# disable predefined suffixes and define your own set of allowed suffixes
 .SUFFIXES:
 .SUFFIXES: .o .f .F .c .f90 .F90

#*************************************************************************
# inference rules (how to compile object files that have no explicit rules)
#  $* = name part of target
#  $@ = full target name
#  $< = dependent name

.f90.o:
	$(CC) -c $(FFLAGS) $<

#**************************************************************************
# Rules for building EXENAME from OBJS and OBJS from your source.

$(EXEOLD):  $(OBJS) $(OBJS2) 
	$(LINK) -o $(EXEOLD) $(OBJS2) 

#************************************************************************
# if you wish to compile a certain object with different flags
# or in some special way, then specify the target & dependency explicitly
# the following line say Timer.o is depended on Timer.f90
#Timer.o: Timer.f90
#	$(CC) -c -O3 Timer.f90

	cp  mod_mathfunc.o mod_polymap.o mod_tpsa.o TPSA 
	cp  AccSimulator.o main.o Input.o Output.o Contrl
	cp  DriftTube.o \
	    ConstFoc.o BeamLineElem.o BeamBunch.o CompDom.o \
	    Distribution.o Appl
	cp  Timer.o Func
	cp  NumConst.o PhysConst.o Data.o Pgrid.o DataStruct 
	cp  ImpactZexe testrun
#***********************************************************************
clean:
	-rm *.o $(EXEOLD) *.mod */*.o
