CXX=g++
CXXFLAGS=-O3 -Wall -Wextra -std=c++11

.PHONY: all clean run plot

all: cache_size

clean:
	rm -f cache_size results.txt results.png

run:
	./cache_size | tee results.txt
	# Results stored into 'results.txt', run 'make plot'.

plot:
	./plot.sh
	# Plot stored into 'results.png'.

cache_size: cache_size.cpp Makefile
	$(CXX) $(CXXFLAGS) cache_size.cpp -o cache_size

