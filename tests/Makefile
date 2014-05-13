BINARIES_DIRECTORY:=bin

SOURCES:=$(wildcard *.cpp)
BINARIES:=$(addprefix $(BINARIES_DIRECTORY)/, $(SOURCES:.cpp=))

#DIRECTIVES:=-O2 -DNDEBUG
DIRECTIVES:=-g

all: $(BINARIES)

$(BINARIES_DIRECTORY)/% : %.cpp $(HEADERS)
	mkdir -p $(BINARIES_DIRECTORY)
	$(CXX) \
			-Wall -pedantic                  \
			-msse2                           \
			-std=c++11 -fopenmp              \
			-I ../../                        \
			$(DIRECTIVES) -o $@ $<

clean:
	$(RM) -r $(BINARIES_DIRECTORY)

.PHONY: all clean

