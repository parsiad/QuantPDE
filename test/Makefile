BIN_DIR:=build

SOURCES:=$(wildcard *.cpp)
BIN:=$(addprefix $(BIN_DIR)/, $(SOURCES:.cpp=))

EIGEN:=-I/usr/include/eigen3 -I/usr/local/include/eigen3
GCC:=/usr/bin/gcc-4.8

OCTAVE_INC_DIR:=/usr/include/octave-3.8.2
OCTAVE_LIB_DIR:=/usr/lib/octave/3.8.2

OPTS:=

DEBUG_OPTS:=-g3
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

all: $(BIN)

$(BIN_DIR)/experimental/hjbqvi: experimental/hjbqvi.cpp
	mkdir -p $(BIN_DIR)
	mkdir -p $(BIN_DIR)/experimental
	$(CXX) \
			-std=c++11 $(OPTS) $(CFLAGS) \
			-I $(OCTAVE_INC_DIR) -I $(OCTAVE_INC_DIR)/octave \
			-L $(OCTAVE_LIB_DIR) -l octinterp -l octave \
			-I ../ $(EIGEN) -o $@ $<

$(BIN_DIR)/% : %.cpp
	mkdir -p $(BIN_DIR)
	mkdir -p $(BIN_DIR)/experimental
	$(CXX) \
			-std=c++11 $(OPTS) $(CFLAGS) \
			-I../ $(EIGEN) -o $@ $<

clean:
	$(RM) -r $(BIN_DIR)

.PHONY: all clean

