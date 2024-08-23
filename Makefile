export CXX = g++
export CXXFLAGS = -g -Wall -fPIC -DNO_ART -O2
export ROOTFLAGS = $(shell root-config --cflags)
export INCLUDE = -I/usr/include/
INCLUDE += -I/uni-mainz.de/homes/losulliv/dune/Genie/build/include/GENIE/
#INCLUDE += -I$(NUSYST) -I$(NUSYST)/build/systematicstools/src/systematicstools
#INCLUDE += -I$(NUSYST)/build/Linux/include/
INCLUDE += -I/usr/include/boost/
INCLUDE += -I/uni-mainz.de/homes/losulliv/dune/
INCLUDE += -I/uni-mainz.de/homes/losulliv/dune/fhiclcpp-simple
#INCLUDE += -I$(CETLIB_INC)
#INCLUDE += -I$(CETLIB_EXCEPT_INC)
INCLUDE += -I/uni-mainz.de/homes/losulliv/dune/duneanaobj-03_03_00/ 

export LDLIBS += -L$(LOG4CPP_LIB) -llog4cpp
LDLIBS += -L$(TBB_LIB) -ltbb
LDLIBS += -L$(LIBXML2_FQ_DIR)/lib -lxml2
LDLIBS += -L$(HDF5_LIB) -lhdf5_cpp
LDLIBS += -L$(PYTHIA6) -lPythia6

DUNEANAOBJ_LIB = /uni-mainz.de/homes/losulliv/dune/duneanaobj-03_03_00/build/duneanaobj/StandardRecord/libduneanaobj_StandardRecord.a
#LDLIBS += -L$(NUSYST)/build/Linux/lib -lsystematicstools_utility -lsystematicstools_interpreters -lsystematicstools_interface -lsystematicstools_systproviders
#LDLIBS += -L$(NUSYST)/build/nusystematics/artless -lnusystematics_systproviders
LDLIBS += -L$(DUNEANAOBJ_LIB) -lduneanaobj_StandardRecord -lduneanaobj_StandardRecordFlat

LDLIBS += -L$(GSL_LIB) -lgsl -lgslcblas
#LDLIBS += -L$(LHAPDF_LIB) -lLHAPDF
LDLIBS += $(shell genie-config --libs)
LDLIBS += $(shell root-config --glibs)
LDLIBS += -lMathMore
LDLIBS += -lGeom -lEGPythia6 -lGenVector
LDLIBS += -Lboost -lboost_program_options
#LDLIBS += -L$(CETLIB_LIB) -L$(CETLIB_EXCEPT_LIB) -lcetlib -lcetlib_except
LDLIBS += -L$(FHICLCPP_LIB) -lfhiclcpp -lfhiclcpp_types

export LIBDIR = $(PWD)/lib
export BINDIR = $(PWD)/bin

SUBDIRS = src

all:
	test -d $(LIBDIR) || mkdir $(LIBDIR)
	test -d $(BINDIR) || mkdir $(BINDIR)
	+make -C src $(MAKECMDGOALS)

clean:
	rm -f $(LIBDIR)/*
	rm -f $(BINDIR)/*
	rm -f $(wildcard AutoDict_*)
	+make -C src clean

test:
	+make -C src $(MAKECMDGOALS)
