CXX = g++
MPICXX = mpicxx

CXXFLAGS = -O2 -g
EXECS = serial mpi mpi_non_blocking

all: $(EXECS)

serial: serial.cpp
	$(CXX) $< -std=c++11 $(CXXFLAGS) -o $@

mpi: mpi.cpp
	$(MPICXX) $< -std=c++11 $(CXXFLAGS) -o $@

mpi_non_blocking: mpi_non_blocking.cpp
	$(MPICXX) $< -std=c++11 $(CXXFLAGS) -o $@

clean:
	rm -f $(EXECS)
