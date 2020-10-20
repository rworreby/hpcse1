CXX = c++
CXXFLAGS = -O3 -Wall -Wextra -Wpedantic

all: main

main: main.cpp Makefile
	$(CXX) -std=c++20 -fopenmp $(CXXFLAGS) main.cpp -o main

run:
	./varym $(N)

plot:
	./plot

runplot: run plot

clean:
	rm -rf main out results.png

.PHONY: all clean run plot runplot
