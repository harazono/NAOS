DEPENDS := .depend.mk
CXXFLAGS := -g -Wall -O3 -DNDEBUG -std=c++11 -fopenmp -lgomp -I /home/harazono/local/include/htslib/
# CXXFLAGS := -g -Wall -O0 -std=c++11 -fopenmp -lgomp
GTESTFLAGS := -lgtest -lpthread
SOURCES := naos.cpp
CXX := $(HOME)/local/bin/g++

.PHONY: depend clean test build_test

all: Naos

Naos: naos.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	-rm *.o

