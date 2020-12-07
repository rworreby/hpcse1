MPICXX = mpicxx
CXXFLAGS = -O2 -g

all: diffusion

diffusion: diffusion.cpp
	$(MPICXX) $< -std=c++11 $(CXXFLAGS) -o $@

run: diffusion
	mpirun -n 1 ./diffusion 1.0 2.0 256

clean:
	rm -f diffusion *.dat

.PHONY: all clean
.SUFFIXES:
