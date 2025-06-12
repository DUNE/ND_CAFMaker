SRCDIR = $(shell pwd)
SOURCES = $(shell find $(SRCDIR) -name '*.cxx')
HEADERS = $(shell find $(SRCDIR) -name '*.h')
OBJS = $(patsubst %.cxx, %.o, $(SOURCES))

INCLUDE += -I$(SRCDIR) -I$(SRPROXY_INC) -I/usr/include/curl
LDLIBS += -L$(SANDRECO_LIB) -lStruct -L/usr/lib64 -lcurl

# enable TMS unless otherwise specified
# (set ENABLE_TMS=0 in the make command line)
ENABLE_TMS := 1
ifneq ($(ENABLE_TMS), 1)
	CXXFLAGS += -DDISABLE_TMS
endif


# # Automatically enable SAND if SANDRECO_INC or SANDRECO_LIB are set
# ifneq ($(origin SANDRECO_INC), undefined)
# 	INCLUDE += -I$(SANDRECO_INC)
# 	SAND_OPT = -DENABLE_SAND
# endif

# ifneq ($(origin SANDRECO_LIB), undefined)
# 	LDLIBS += -L$(SANDRECO_LIB) -lStruct
# 	SAND_OPT = -DENABLE_SAND
# endif


LIBNAME = ND_CAFMaker
LIBS = $(LIBDIR)/lib$(LIBNAME).so
BINS = $(BINDIR)/makeCAF

# Option to enable struct.h dictionary (default: off)
ENABLE_SAND_DICT ?= 0

# Add variables for struct dictionary and library
DEPSDIR = $(SRCDIR)/deps
STRUCT_HDR = $(DEPSDIR)/struct.h
STRUCT_LINKDEF = $(DEPSDIR)/StructLinkDef.h
STRUCT_DICT_SRC = $(LIBDIR)/StructDict.cxx
STRUCT_DICT_OBJ = $(LIBDIR)/StructDict.o
STRUCT_LIB = $(LIBDIR)/libStruct.so
STRUCT_PCM = $(LIBDIR)/libStruct_rdict.pcm

ifeq ($(ENABLE_SAND_DICT),1)
$(shell wget https://raw.githubusercontent.com/DUNE/sandreco/refs/heads/develop/include/struct.h -P $(DEPSDIR))
$(shell wget https://raw.githubusercontent.com/DUNE/sandreco/refs/heads/develop/include/StructLinkDef.h -P $(DEPSDIR))
INCLUDE += -I$(DEPSDIR)
SAND_OPT = -DENABLE_SAND
LDLIBS_EXEC += -L$(LIBDIR) -lStruct
DICT_TARGETS = $(STRUCT_LIB)
else
DICT_TARGETS =
LDLIBS_EXEC =
SAND_OPT =
endif

RPATH_FLAG = -Wl,-rpath,$(DEPSDIR)

# build the library and the binaries
all: $(BINS)

# testing...
test: $(BINDIR)/testHDF

# make an object for any .cxx requested
%.o : %.cxx $(if $(wildcard %.h),%.h,)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCLUDE) -fPIC $(SAND_OPT) -o $@ -c $< #compile

# should *really* make a .d rule for the objects
# to get their dependencies on the headers right, but ... not now
$(LIBDIR)/lib$(LIBNAME).so: $(OBJS) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -shared $(RPATH_FLAG) -o $@ $(OBJS)

$(BINDIR)/makeCAF: $(LIBS) makeCAF.C $(DICT_TARGETS)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCLUDE) -L$(LIBDIR) -l$(LIBNAME) $(LDLIBS) $(LDLIBS_EXEC) $(RPATH_FLAG) -o $@ makeCAF.C

$(BINDIR)/testHDF: $(LIBS) testHDF.C $(DICT_TARGETS)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCLUDE) -L$(LIBDIR) -l$(LIBNAME) $(LDLIBS) $(LDLIBS_EXEC) $(RPATH_FLAG) -o $@ testHDF.C

ifeq ($(ENABLE_SAND_DICT),1)

$(STRUCT_DICT_SRC): $(STRUCT_HDR) $(STRUCT_LINKDEF)
	$(ROOTSYS)/bin/rootcling -f $(STRUCT_DICT_SRC) -s $(STRUCT_LIB) -rml libStruct.so -c $(INCLUDE) $^

$(STRUCT_DICT_OBJ): $(STRUCT_DICT_SRC)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCLUDE) -fPIC -o $@ -c $<

# Rule to build the shared library for struct
$(STRUCT_LIB): $(STRUCT_DICT_OBJ)
	$(CXX) -shared -o $@ $^
endif

clean:
	rm -f $(OBJS)
	rm -f $(LIBS)
	rm -f $(BINS)
	rm -f $(wildcard AutoDict_*)
	rm -f $(LIBDIR)/StructDict.cxx $(LIBDIR)/StructDict.o $(LIBDIR)/libStruct.so $(LIBDIR)/libStruct_rdict.pcm
