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
CXXFLAGS = -g -DGCC4 -DPLUGGABLE_SCORING $(OSFLAGS) $(ZLIB_INCL) -DREAD_GZIP -DHAVE_ZLIB -I ../../gzstream $(PWIZ_MZML_FLAGS) 
#CXXFLAGS = -O2 -DGCC -DPLUGGABLE_SCORING -DX_P3  $(ZLIB_INCL)
ifeq (${OS},Windows_NT)
# MinGW
LDFLAGS = -static -static-libstdc++ 
else
LDFLAGS = -lpthread
endif

# optional mac os x library location
LDFLAGS += -L/usr/lib
ifeq (exists, $(shell [ -d /sw/lib ] ) && echo exists )
LDFLAGS += -L/sw/lib
endif

LDFLAGS += -lm $(EXPAT_LIB) $(ZLIB_LIB)

# building multinode version (X!!Tandem)?
ifneq ("$(XBANGBANG)","")
EXECUTABLE = $(BUILD_DIR)/bbtandem
XVARIANT = _xbangbang_
LINKCC = mpicxx
CXXFLAGS += -DXBANGBANG -DHAVE_MULTINODE_TANDEM
# [multinode] MPI: look for 64-bit mpich2 libraries and headers
ifneq "$(wildcard /usr/lib64/mpich2/lib/libmpichcxx.a )" ""
CXXFLAGS += -I /usr/include/mpich2-x86_64
LDFLAGS += /usr/lib64/mpich2/lib/libmpichcxx.a /usr/lib64/mpich2/lib/libmpich.a
else
# [multinode] MPI: look for other mpich2 libraries and headers
ifneq "$(wildcard /usr/local/mpich2/include )" ""
CXXFLAGS += -I /usr/local/mpich2/include -I /usr/local/mpich2
LDFLAGS += /usr/local/mpich2/lib/libmpichcxx.a /usr/local/mpich2/lib/libmpich.a
else # not using mpich2. look for open mpi
# [multinode] MPI: look for open mpi headers as in StarCluster ubuntu;
# if found also use custom compiler
ifneq "$(wildcard /usr/lib/openmpi/include/mpi.h )" ""
CXXFLAGS += -I /usr/lib/openmpi/include
CXX = mpicxx
endif
# [multinode] MPI: look for 64-bit open mpi headers as in StarCluster
# centos; if found also use custom compiler
ifneq "$(wildcard /usr/lib64/openmpi/1.4-gcc/include/mpi.h )" ""
CXXFLAGS += -I /usr/lib64/openmpi/1.4-gcc/include
CXX = mpicxx
endif
# [multinode] MPI: look for standard open mpi headers as in
# StarCluster centos; if found also use custom compiler
ifneq "$(wildcard /usr/lib/openmpi/1.4-gcc/include/mpi.h )" ""
CXXFLAGS += -I /usr/lib/openmpi/1.4-gcc/include
CXX = mpicxx
endif
endif # end open search (vs mpich2)
endif
endif # end multinode build (X!!Tandem)

# building our GPU-enabled X!tandem build, tandem-g?
ifneq ("$(XTGPU)","")
XVARIANT = _xgpu_
EXECUTABLE = $(BUILD_DIR)/tandem-g
NVCC_WARNINGS = --compiler-options -Wall
# -Wall disabled for now, due to lots of messy 64-bit warnings and
# -unused variable in original code; otherwise use:
# NVCC_WARNINGS = --compiler-options -Wall --compiler-options -Werror

# TODO: autodetect this
NVCC_TARGET = --machine 64

NVCC_DEBUG = -g
NVCC_PROFILE = -pg

# (Tip: add -keep to NVCC in order to preserve intermediate files
# (.cu.cpp, .ptx, etc) during development)

NVCC=nvcc $(NVCC_DEBUG) $(NVCC_PROFILE)  $(NVCC_WARNINGS) $(NVCC_TARGET) $(CCOPT) $(INCLUDE)
NVCC_LIBS =-L/usr/local/cuda/lib -lcuda -lcublas -lcudart

# attempt to get rid of "rpath" issue for debugging purposes (os x only)
ifeq ($(ARCH_FAMILY),darwin)
NVCC_DEBUG += -Xlinker -rpath -Xlinker /usr/local/cuda/lib 
endif

LINKCC = $(CXX)
endif
# (end of GPGPU build)


ifeq ("$(XTGPU)$(XBANGBANG)","")
# normal TPP X!Tandem build
LINKCC = $(CXX)
endif

# certain specific files require the GPGPU compiler for GPU builds and
# use the standard compiler otherwise
ifeq "$(XVARIANT)" "_xgpu_"
XVARIANT_CC=$(NVCC)
LINKCC=$(NVCC)
else
XVARIANT_CC=$(CXX)
endif


.SUFFIXES:	.o .cpp

OBJCLASS = $(ARCH)$(XVARIANT)

SRCS := $(wildcard *.cpp)
# (make sure GPGPU .cpp files don't get added to the source file list
# for non GPGPU build variants)
ifneq "$(XVARIANT)" "_xgpu_"
SRCS_FILTERED := $(filter-out mscore_kgpu.cpp, $(SRCS))
SRCS := $(SRCS_FILTERED)
endif

OBJS := $(patsubst %.cpp,$(OBJCLASS)%.o,$(SRCS))

# for GPGPU build variant, add additional GPGPU files, which are always
# built with the CUDA compiler
ifeq "$(XVARIANT)" "_xgpu_"
SRCS += $(wildcard *.cu)
OBJS += $(patsubst %.cu,$(OBJCLASS)%.o,$(wildcard *.cu))
endif



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
	$(LINKCC) $(CXXFLAGS) -o $(EXECUTABLE) $(OBJS) $(PWIZ_LIBS)  $(LDFLAGS) $(GZSTREAMLIB) $(ZLIB_LIB) $(USER_OBJS) $(NVCC_LIBS)

# specify rules for creating the dependency files, which depend on the cpp files

# (this file is only built in the GPGPU variant, using the CUDA compiler, but gets included in the dependency wildcard for all builds so we might generated dependencies with the standard compiler)
$(OBJCLASS)mscore_kgpu.d: mscore_kgpu.cpp
	$(XVARIANT_CC) -M $(CXXFLAGS) $< > $@
	$(XVARIANT_CC) -M $(CXXFLAGS) $< | sed s/\\.o/.d/ > $@

# these files are only built in the GPGPU variant, using the CUDA compiler
$(OBJCLASS)%.d: %.cu
	$(NVCC) -M $(CXXFLAGS) $< > $@
	$(NVCC) -M $(CXXFLAGS) $< | sed s/\\.o/.d/ > $@

# this file is only built in the GPGPU variant, using the CUDA compiler
#$(OBJCLASS)mscore_kgpu_thrust.d: mscore_kgpu_thrust.cu
#	$(NVCC) -M $(CXXFLAGS) $< > $@
#	$(NVCC) -M $(CXXFLAGS) $< | sed s/\\.o/.d/ > $@

# all other files use the specified compiler for all builds
$(OBJCLASS)%.d: %.cpp
	$(XVARIANT_CC) -M $(CXXFLAGS) $< > $@
	$(XVARIANT_CC) -M $(CXXFLAGS) $< | sed s/\\.o/.d/ > $@




# specify rules for the object files

$(OBJCLASS)mscore_kgpu.o: mscore_kgpu.cpp
	$(XVARIANT_CC) $(CXXFLAGS) -c -o $@ $<

# this file is only built in the GPGPU build, using the CUDA compiler
$(OBJCLASS)mscore_kgpu_thrust.o: mscore_kgpu_thrust.cu
	$(NVCC) $(CXXFLAGS) -c -o $@ $<


$(OBJCLASS)%.o: %.cpp
	$(XVARIANT_CC) $(CXXFLAGS) -c -o $@ $<

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
