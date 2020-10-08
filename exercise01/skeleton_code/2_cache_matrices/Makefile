CXX = g++
CXXFLAGS = -O3 -std=c++11 -Wall -Wextra

targets = matrix_vector transpose  matrix_matrix

all: $(targets)

% : %.cpp
	$(CXX) $< $(CXXFLAGS) -o $@

clean:
	rm -f *.o *~ $(targets) *.txt

.PHONY: all
.PHONY: clean
