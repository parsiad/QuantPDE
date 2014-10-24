BINARIES_DIRECTORY:=build

SOURCES:=$(wildcard *.cpp)
BINARIES:=$(addprefix $(BINARIES_DIRECTORY)/, $(SOURCES:.cpp=))

EIGEN:=-I/usr/local/include/eigen3 -I/usr/include/eigen3
GCC:=/usr/bin/gcc-4.8

ifdef DEBUG
DIRECTIVES:=-static -g3 -DEIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
else
DIRECTIVES:=-O2 -DNDEBUG
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

