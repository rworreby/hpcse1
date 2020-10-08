CXX = c++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -Wpedantic

targets = matrix_vector transpose  matrix_matrix

.PHONY: all
all: $(targets)

% : %.cpp
	$(CXX) $< $(CXXFLAGS) -o $@

.PHONY: clean
clean:
	rm -f *.o *~ $(targets) *.txt
