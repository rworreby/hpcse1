MPICXX = mpicxx
CXXFLAGS = -O2 -g

all: diffusion

diffusion: diffusion.cpp
	$(MPICXX) $< -std=c++11 $(CXXFLAGS) -o $@

clean:
	rm -f diffusion *.dat

.PHONY: all clean
.SUFFIXES:
