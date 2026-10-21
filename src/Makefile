SRCDIR = $(shell pwd)
SOURCES = $(shell find $(SRCDIR) -name '*.cxx')
HEADERS = $(shell find $(SRCDIR) -name '*.h')
OBJS = $(patsubst %.cxx, %.o, $(SOURCES))

INCLUDE += -I$(SRCDIR) -I$(SRPROXY_INC) -I/usr/include/curl
LDLIBS += -L/usr/lib64 -lcurl

ENABLE_TMS := 1
ifneq ($(ENABLE_TMS),1)
  CXXFLAGS += -DDISABLE_TMS
endif

LIBNAME = ND_CAFMaker
LIBS    = $(LIBDIR)/lib$(LIBNAME).so
BINS    = $(BINDIR)/makeCAF

ENABLE_SAND_DICT ?= 0

SANDRECO_INC := $(SRCDIR)/../sandreco/sandreco/include

STRUCT_HDR        := $(SANDRECO_INC)/struct.h
STRUCT_LINKDEF    := $(SANDRECO_INC)/StructLinkDef.h
STRUCT_DICT_SRC   := $(LIBDIR)/StructDict.cxx
STRUCT_DICT_OBJ   := $(LIBDIR)/StructDict.o
STRUCT_LIB        := $(LIBDIR)/libStruct.so
STRUCT_PCM        := $(LIBDIR)/libStruct_rdict.pcm

ifeq ($(ENABLE_SAND_DICT),1)
  $(shell mkdir -p $(LIBDIR))
  INCLUDE      += -I$(SANDRECO_INC)
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

%.o : %.cxx $(if $(wildcard %.h),%.h,)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) $(INCLUDE) -fPIC $(SAND_OPT) -o $@ -c $<

$(LIBDIR)/lib$(LIBNAME).so: $(OBJS)
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

clean:
	rm -f $(OBJS)
	rm -f $(LIBS)
	rm -f $(BINS)
	rm -f $(STRUCT_DICT_SRC) $(STRUCT_DICT_OBJ)
	rm -f $(STRUCT_LIB) $(STRUCT_PCM)