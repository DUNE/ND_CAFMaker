SRCDIR = $(shell pwd)
SOURCES = $(shell find $(SRCDIR) -name '*.cxx')
HEADERS = $(shell find $(SRCDIR) -name '*.h')
OBJS = $(patsubst %.cxx, %.o, $(SOURCES))

INCLUDE += -I$(SRCDIR)


# Automatically enable SAND if SANDRECO_INC or SANDRECO_LIB are set
ifneq (x$(SANDRECO_INC), x)
	INCLUDE += -I$(SANDRECO_INC)
	SAND_OPT = -DENABLE_SAND
endif

ifneq (x$(SANDRECO_LIB), x)
	LDLIBS += -L$(SANDRECO_LIB) -lStruct
	SAND_OPT = -DENABLE_SAND
endif


LIBNAME = ND_CAFMaker
LIBS = $(LIBDIR)/lib$(LIBNAME).so
BINS = $(BINDIR)/makeCAF

# build the library and the binaries
all: $(BINS)

# testing...
test: $(BINDIR)/testHDF

# make an object for any .cxx requested
%.o : %.cxx %.h
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) -fPIC $(SAND_OPT) -o $@ -c $< #compile

# should *really* make a .d rule for the objects
# to get their dependencies on the headers right, but ... not now
$(LIBDIR)/lib$(LIBNAME).so: $(OBJS) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -shared -o $@ $(OBJS)

$(BINDIR)/makeCAF: $(LIBS) makeCAF.C
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) $(LDLIBS) -L$(LIBDIR) -l$(LIBNAME) -o $@ makeCAF.C

$(BINDIR)/testHDF: $(LIBS) testHDF.C
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) $(LDLIBS) -L$(LIBDIR) -l$(LIBNAME) -o $@ testHDF.C


clean:
	rm -f $(OBJS)
	rm -f $(LIBS)
	rm -f $(BINS)
	rm -f $(wildcard AutoDict_*)
