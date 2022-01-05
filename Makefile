export CXX = g++
export CXXFLAGS = -g -Wall -fPIC -DNO_ART -O2
export ROOTFLAGS = `root-config --cflags --glibs`
export INCLUDE = -I$(HDF5_INC)
INCLUDE += -I$(GENIE_INC)/GENIE
#INCLUDE += -I$(NUSYST) -I$(NUSYST)/build/systematicstools/src/systematicstools
#INCLUDE += -I$(NUSYST)/build/Linux/include/
INCLUDE += -I$(DUNEANAOBJ_INC)

export LDLIBS += -L$(LOG4CPP_LIB) -llog4cpp
LDLIBS += -L/usr/lib64 -lxml2
LDLIBS += -L$(HDF5_LIB) -lhdf5_cpp
LDLIBS += -L$(PYTHIA6) -lPythia6
LDLIBS += -L$(ROOTSYS)/lib -lGeom -lEGPythia6
LDLIBS += -L$(GENIE)/lib \
                                        -lGAlgorithm \
                                        -lGBaryonResonance \
                                        -lGBase \
                                        -lGBodekYang \
                                        -lGCharm \
                                        -lGCoh \
                                        -lGCrossSections \
                                        -lGDecay \
                                        -lGDfrc \
                                        -lGDIS \
                                        -lGElas \
                                        -lGElFF \
                                        -lGEVGCore \
                                        -lGEVGDrivers \
                                        -lGEVGModules \
                                        -lGFluxDrivers \
                                        -lGFragmentation \
                                        -lGGeo \
                                        -lGGiBUU \
                                        -lGHadronTransp \
                                        -lGHEP \
                                        -lGInteraction \
                                        -lGLlewellynSmith \
                                        -lGMEC \
                                        -lGReinSehgal \
                                        -lGSingleKaon \
                                        -lGMessenger \
                                        -lGMuELoss \
                                        -lGNtuple \
                                        -lGNuclear \
                                        -lGNuE \
                                        -lGNuGamma \
                                        -lGNumerical \
                                        -lGPDF \
                                        -lGPDG \
                                        -lGQEL \
                                        -lGQPM \
                                        -lGRegistry \
                                        -lGRES \
                                        -lGUtils \
                                        -lGReWeight

#LDLIBS += -L$(NUSYST)/build/Linux/lib -lsystematicstools_utility -lsystematicstools_interpreters -lsystematicstools_interface -lsystematicstools_systproviders
#LDLIBS += -L$(NUSYST)/build/nusystematics/artless -lnusystematics_systproviders
LDLIBS += -L$(DUNEANAOBJ_LIB) -lduneanaobj_StandardRecord

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