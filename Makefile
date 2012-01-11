#makefile for c++ programs
#some tweaks by Brian Pratt for TPP usage, esp. multiplatform builds
#even more tweaks for MapReduce and MPI (X!!Tandem) implementations

# make it possible to build mingw and cygwin on same machine
ifeq ($(BUILD_DIR),)
# no build dir specified, generic usage
BUILD_DIR=../bin
else
# for TPP install we want these files as well
SUPPORTFILES= $(BUILD_DIR)/tandem_params.xml $(BUILD_DIR)/isb_default_input_kscore.xml $(BUILD_DIR)/isb_default_input_native.xml $(BUILD_DIR)/taxonomy.xml
# for pwiz mzML read lib
include $(SRC_ROOT)Makefile.incl
PWIZ_MZML_FLAGS= -DHAVE_PWIZ_MZML_LIB $(PWIZ_INCL) -I $(SRC_ROOT)Parsers/ramp 
PWIZ_LIBS= $(OBJS_PWIZ) $(BOOST_LIBS)  
VPATH= $(VPATH_PWIZ)
endif

EXECUTABLE = $(BUILD_DIR)/tandem
GZSTREAMLIB = $(BUILD_DIR)/libgzstream.a

ESCAPED_TANDEM_PARAMETERS_INSTALL_DIR=$(subst \,\\,$(TANDEM_PARAMETERS_INSTALL_DIR))


#EXECUTABLE = $(BUILD_DIR)/p3.exe


#CXXFLAGS denotes flags for the C++ compiler

CXX = g++
CXXFLAGS = -O2 -DGCC4 -DPLUGGABLE_SCORING $(OSFLAGS) $(ZLIB_INCL) -DREAD_GZIP -DHAVE_ZLIB -DHAVE_MULTINODE_TANDEM -I ../../gzstream $(PWIZ_MZML_FLAGS) 
#CXXFLAGS = -O2 -DGCC -DPLUGGABLE_SCORING -DX_P3  $(ZLIB_INCL)
ifeq (${OS},Windows_NT)
# MinGW
LDFLAGS = -static -static-libstdc++ 
else
LDFLAGS = -lpthread
endif
LDFLAGS += -L/usr/lib -L/sw/lib -lm $(EXPAT_LIB) $(ZLIB_LIB)

# building for MPI?
ifneq ("$(XBANGBANG)","")
XVARIANT = _xbangbang_
LINKCC = mpicxx
CXXFLAGS += -DXBANGBANG
# is it mpich2?
ifneq "$(wildcard /usr/lib64/mpich2/lib/libmpichcxx.a )" ""
CXXFLAGS += -I /usr/include/mpich2-x86_64
LDFLAGS += /usr/lib64/mpich2/lib/libmpichcxx.a /usr/lib64/mpich2/lib/libmpich.a
else
ifneq "$(wildcard /usr/local/mpich2/include )" ""
CXXFLAGS += -I /usr/local/mpich2/include -I /usr/local/mpich2
LDFLAGS += /usr/local/mpich2/lib/libmpichcxx.a /usr/local/mpich2/lib/libmpich.a
else
# is it open mpi?
# as in StarCluster ubuntu
ifneq "$(wildcard /usr/lib/openmpi/include/mpi.h )" ""
CXXFLAGS += -I /usr/lib/openmpi/include
CXX = mpicxx
endif
# as in StarCluster centos
ifneq "$(wildcard /usr/lib64/openmpi/1.4-gcc/include/mpi.h )" ""
CXXFLAGS += -I /usr/lib64/openmpi/1.4-gcc/include
CXX = mpicxx
endif
# as in StarCluster centos
ifneq "$(wildcard /usr/lib/openmpi/1.4-gcc/include/mpi.h )" ""
CXXFLAGS += -I /usr/lib/openmpi/1.4-gcc/include
CXX = mpicxx
endif
endif
endif
EXECUTABLE = $(BUILD_DIR)/bbtandem

else
LINKCC = $(CXX)
endif



.SUFFIXES:	.o .cpp

OBJCLASS = $(ARCH)$(XVARIANT)

SRCS := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,$(OBJCLASS)%.o,$(wildcard *.cpp))
DEPS := $(patsubst %.o,%.d,$(OBJS))

all: $(EXECUTABLE) $(SUPPORTFILES)

# define some support files for the win32 TPP installer
$(BUILD_DIR)/tandem_params.xml: ../bin/tandem_params.xml
	cp $< $@
	echo escaped: ${ESCAPED_TANDEM_PARAMETERS_INSTALL_DIR}
	sed 's?_DEFAULT_INPUT_LOCATION_/?${ESCAPED_TANDEM_PARAMETERS_INSTALL_DIR}?g' $@ > ${@}.updated
	mv ${@}.updated $@

$(BUILD_DIR)/isb_default_input_kscore.xml: ../bin/isb_default_input_kscore.xml
	cp $< $@

$(BUILD_DIR)/isb_default_input_native.xml: ../bin/isb_default_input_native.xml
	cp $< $@

$(BUILD_DIR)/taxonomy.xml: ../bin/taxonomy.xml
	cp $< $@

#define the components of the program, and how to link them
#these components are defined as dependencies; that is they must be up-to-date before the code is linked

$(EXECUTABLE): $(DEPS) $(OBJS) $(GZSTREAMLIB) $(PWIZ_LIBS) $(USER_OBJS)
	$(LINKCC) $(CXXFLAGS) -o $(EXECUTABLE) $(OBJS)  $(PWIZ_LIBS)  $(LDFLAGS) $(GZSTREAMLIB) $(ZLIB_LIB) $(USER_OBJS)

#specify the dep files depend on the cpp files

$(OBJCLASS)%.d: %.cpp
	$(CXX) -M $(CXXFLAGS) $< > $@
	$(CXX) -M $(CXXFLAGS) $< | sed s/\\.o/.d/ > $@

$(OBJCLASS)%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(EXECUTABLE) $(DEPS) $(SUPPORTFILES)

explain:
	@echo "The following info represents the program:"
	@echo "Final exec name: $(EXECUTABLE)"
	@echo "Source files:       $(SRCS)"
	@echo "Object files:       $(OBJS)"
	@echo "Dep files:          $(DEPS)"

depend: $(DEPS)
	@echo "Deps are now up-to-date."
 
-include $(DEPS)
