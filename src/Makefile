SRCDIR = $(shell pwd)
SOURCES = $(shell find $(SRCDIR) -name '*.cxx')
HEADERS = $(shell find $(SRCDIR) -name '*.h')
OBJS = $(patsubst %.cxx, %.o, $(SOURCES))

INCLUDE += -I$(SRCDIR)

LIBNAME = ND_CAFMaker
LIBS = $(LIBDIR)/lib$(LIBNAME).so
BINS = $(BINDIR)/makeCAF

# build the library and the binaries
all: $(BINS)

# make an object for any .cxx requested
%.o : %.cxx
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) -fPIC -o $@ -c $^ #compile

# should *really* make a .d rule for the objects
# to get their dependencies on the headers right, but ... not now
$(LIBDIR)/lib$(LIBNAME).so: $(OBJS) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -shared -o $@ $(OBJS)

$(BINDIR)/makeCAF: $(LIBS) makeCAF.C
	@echo "making makeCAF"
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) $(LDLIBS) -L$(LIBDIR) -l$(LIBNAME) -o $@ makeCAF.C


clean:
	rm -f $(OBJS)
	rm -f $(LIBS)
	rm -f $(BINS)
	rm -f $(wildcard AutoDict_*)