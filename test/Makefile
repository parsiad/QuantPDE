BINARIES_DIRECTORY:=build

SOURCES:=$(wildcard *.cpp)
BINARIES:=$(addprefix $(BINARIES_DIRECTORY)/, $(SOURCES:.cpp=))

EIGEN:=-I/usr/local/include/eigen3 -I/usr/include/eigen3
GCC:=/usr/bin/gcc-4.8

DEBUG_DIRECTIVES:=-static -g3
PRODUCTION_DIRECTIVES:=-O2 -DNDEBUG

ifdef PROFILE
DIRECTIVES:=$(DEBUG_DIRECTIVES) $(PRODUCTION_DIRECTIVES)
ifdef DEBUG
DIRECTIVES:=$(DEBUG_DIRECTIVES)
else
DIRECTIVES:=$(PRODUCTION_DIRECTIVES)
endif
endif

ifdef CUDA
CXX:=nvcc
OPTIONS:=-x cu --compiler-bindir $(GCC)
else
ifdef CLANG
CXX:=clang++
endif
OPTIONS:=-Wall -pedantic -msse2
endif

all: $(BINARIES)

$(BINARIES_DIRECTORY)/% : %.cpp
	mkdir -p $(BINARIES_DIRECTORY)
	$(CXX) \
			$(OPTIONS) $(CFLAGS) -std=c++11 \
			-I../ $(EIGEN) $(DIRECTIVES) -o $@ $<

clean:
	$(RM) -r $(BINARIES_DIRECTORY)

.PHONY: all clean

