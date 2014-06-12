BINARIES_DIRECTORY:=build

SOURCES:=$(wildcard *.cpp)
BINARIES:=$(addprefix $(BINARIES_DIRECTORY)/, $(SOURCES:.cpp=))

EIGEN:=/usr/local/include/eigen3

ifdef DEBUG
DIRECTIVES:=-g
else
DIRECTIVES:=-O2 -DNDEBUG
endif

all: $(BINARIES)

$(BINARIES_DIRECTORY)/% : %.cpp $(HEADERS)
	mkdir -p $(BINARIES_DIRECTORY)
	$(CXX) \
			-Wall -pedantic                  \
			-msse2                           \
			-std=c++11                       \
			-I ../                           \
			-I $(EIGEN) $(DIRECTIVES) -o $@ $<

clean:
	$(RM) -r $(BINARIES_DIRECTORY)

.PHONY: all clean

