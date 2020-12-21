CXX = g++
CXXFLAGS = -Wall -O3 -g -std=c++11 -DNDEBUG

all: diffusionADI

diffusionADI: diffusionADI.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f diffusionADI *.o

cleandat:
	rm -f *.dat *.png

cleanall: clean cleandat

.PHONY: clean cleanall cleandat all
.SUFFIXES:
