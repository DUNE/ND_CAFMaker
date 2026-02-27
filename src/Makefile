SRCDIR = $(shell pwd)
SOURCES = $(shell find $(SRCDIR) -name '*.cxx')
HEADERS = $(shell find $(SRCDIR) -name '*.h')
OBJS = $(patsubst %.cxx, %.o, $(SOURCES))

INCLUDE += -I$(SRCDIR) -I$(SRPROXY_INC) -I/usr/include/curl
LDLIBS += -L$(SANDRECO_LIB) -lStruct -L/usr/lib64 -lcurl

ENABLE_TMS := 1
ifneq ($(ENABLE_TMS),1)
  CXXFLAGS += -DDISABLE_TMS
endif

LIBNAME = ND_CAFMaker
LIBS    = $(LIBDIR)/lib$(LIBNAME).so
BINS    = $(BINDIR)/makeCAF

ENABLE_SAND_DICT ?= 0

DEPSDIR           = $(SRCDIR)/deps
STRUCT_HDR        = $(DEPSDIR)/struct.h
STRUCT_LINKDEF    = $(DEPSDIR)/StructLinkDef.h
STRUCT_DICT_SRC   = $(LIBDIR)/StructDict.cxx
STRUCT_DICT_OBJ   = $(LIBDIR)/StructDict.o
STRUCT_LIB        = $(LIBDIR)/libStruct.so
STRUCT_PCM        = $(LIBDIR)/libStruct_rdict.pcm


SANDRECO_TAG  := v02_00_00
SANDRECO_BASE := https://raw.githubusercontent.com/DUNE/sandreco/$(SANDRECO_TAG)/include

ifeq ($(ENABLE_SAND_DICT),1)

  $(shell mkdir -p $(DEPSDIR))
  $(shell wget -N $(SANDRECO_BASE)/struct.h -P $(DEPSDIR))
  $(shell wget -N $(SANDRECO_BASE)/StructLinkDef.h -P $(DEPSDIR))

  INCLUDE      += -I$(DEPSDIR)
  SAND_OPT      = -DENABLE_SAND
  LDLIBS_EXEC  += -L$(LIBDIR) -lStruct
  DICT_TARGETS  = $(STRUCT_LIB)

else

  DICT_TARGETS  =
  LDLIBS_EXEC   =
  SAND_OPT      =

endif

RPATH_FLAG = -Wl,-rpath,$(LIBDIR)

all: $(BINS)

test: $(BINDIR)/testHDF

%.o : %.cxx $(if $(wildcard %.h),%.h,)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCLUDE) -fPIC $(SAND_OPT) -o $@ -c $<


$(LIBDIR)/lib$(LIBNAME).so: $(OBJS) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -shared $(RPATH_FLAG) -o $@ $(OBJS)

ifeq ($(ENABLE_SAND_DICT),1)

  $(STRUCT_DICT_SRC): $(STRUCT_HDR) $(STRUCT_LINKDEF)
	$(ROOTSYS)/bin/rootcling -f $@ -s $(STRUCT_LIB) -rml libStruct.so -c $(INCLUDE) $^

  $(STRUCT_DICT_OBJ): $(STRUCT_DICT_SRC)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCLUDE) -fPIC -o $@ -c $<

  $(STRUCT_LIB): $(STRUCT_DICT_OBJ)
	$(CXX) -shared -o $@ $^

endif


$(BINDIR)/makeCAF: $(LIBS) makeCAF.C $(DICT_TARGETS)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCLUDE) \
	  -L$(LIBDIR) -l$(LIBNAME) \
	  $(LDLIBS) $(LDLIBS_EXEC) \
	  $(RPATH_FLAG) \
	  -o $@ makeCAF.C

$(BINDIR)/testHDF: $(LIBS) testHDF.C $(DICT_TARGETS)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCLUDE) \
	  -L$(LIBDIR) -l$(LIBNAME) \
	  $(LDLIBS) $(LDLIBS_EXEC) \
	  $(RPATH_FLAG) \
	  -o $@ testHDF.C

# should *really* make a .d rule for the objects 
# to get their dependencies on the headers right, but ... not now 


clean:
	rm -f $(OBJS)
	rm -f $(LIBS)
	rm -f $(BINS)
	rm -f $(wildcard AutoDict_*)
	rm -f $(STRUCT_DICT_SRC) $(STRUCT_DICT_OBJ)
	rm -f $(STRUCT_LIB) $(STRUCT_PCM)
	rm -f $(STRUCT_HDR) $(STRUCT_LINKDEF)