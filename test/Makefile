BINARIES_DIRECTORY:=build

SOURCES:=$(wildcard *.cpp)
BINARIES:=$(addprefix $(BINARIES_DIRECTORY)/, $(SOURCES:.cpp=))

EIGEN:=-I/usr/include/eigen3 -I/usr/local/include/eigen3
GCC:=/usr/bin/gcc-4.8

OPTS:=

DEBUG_OPTS:=-static -g3
OPTIMIZE_OPTS:=-O2 -DNDEBUG

ifdef PROFILE
OPTS:=$(OPTS) $(DEBUG_OPTS) $(OPTIMIZE_OPTS)
else
ifdef DEBUG
OPTS:=$(OPTS) $(DEBUG_OPTS)
else
OPTS:=$(OPTS) $(OPTIMIZE_OPTS)
endif
endif

ifdef CUDA
CXX:=nvcc
OPTS:=$(OPTS) -x cu --compiler-bindir $(GCC)
else
ifdef OPENCL
OPTS:=$(OPTS) -DVIENNACL_WITH_OPENCL -lOpenCL
endif
ifdef CLANG
CXX:=clang++
endif
OPTS:=$(OPTS) -Wall -pedantic
endif

all: $(BINARIES)

$(BINARIES_DIRECTORY)/% : %.cpp
	mkdir -p $(BINARIES_DIRECTORY)
	$(CXX) \
			-std=c++11 $(OPTS) $(CFLAGS) \
			-I../ $(EIGEN) -o $@ $<

clean:
	$(RM) -r $(BINARIES_DIRECTORY)

.PHONY: all clean

