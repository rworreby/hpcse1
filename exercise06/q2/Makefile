MPICXX = mpic++
CXXFLAGS = -O2 -g -Wall -Wextra -DVERIFY -fopenmp
LDFLAGS = -lopenblas
BINARIES = cannonMPI cannonHybrid
MPIRUN = mpirun

all: $(BINARIES)

cannonMPI: cannon.cpp
	$(MPICXX) $< -std=c++11 $(CXXFLAGS) $(LDFLAGS) -o $@

cannonHybrid: cannon.cpp
	$(MPICXX) $< -std=c++11 $(CXXFLAGS) $(LDFLAGS) -o $@ -DHYBRID

clean:
	$(RM) $(BINARIES)

test: $(BINARIES)
	@echo "Tests for Cannon (MPI only)"
	$(MPIRUN) -n 1 ./cannonMPI 1024
	$(MPIRUN) -n 4 ./cannonMPI 1024
	$(MPIRUN) -n 16 ./cannonMPI 1024
	@echo "All MPI-only tests passed, you may submit this part of the solution now."
	@echo
	@echo "Tests for Cannon (MPI+OpenMP)"
	OMP_NUM_THREADS=1 $(MPIRUN) -n 1 ./cannonHybrid 512
	OMP_NUM_THREADS=4 $(MPIRUN) -n 1 ./cannonHybrid 512
	OMP_NUM_THREADS=4 $(MPIRUN) -n 4 ./cannonHybrid 512
	OMP_NUM_THREADS=16 $(MPIRUN) -n 4 ./cannonHybrid 512
	OMP_NUM_THREADS=16 $(MPIRUN) -n 16 ./cannonHybrid 512
	@echo "All MPI+OpenMP tests passed, you may submit this part of the solution now."

.PHONY: all clean test
.SUFFIXES:
