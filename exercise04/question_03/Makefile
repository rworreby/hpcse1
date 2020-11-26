# Copyright 2020 ETH Zurich. All Rights Reserved.

# TODO a): Set the compiler
CXX= mpic++
CXXFLAGS+= -std=c++11 -Wall -Wpedantic -Wextra -O3

main: main.o
	$(CXX) $(CXXFLAGS) main.o -o main

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf *.o main
