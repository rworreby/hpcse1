# File       : Makefile
# Description: Compile targets
# Copyright 2020 ETH Zurich. All Rights Reserved.
CXX ?= g++
###############################################################################
# TODO: Fill in the compiler flags.  Make sure you compile SSE code
CXXFLAGS =
###############################################################################

# debug=true --> use assertions for debugging
debug ?= false

.PHONY: clean measurement

ifeq ("$(debug)", "false")
###############################################################################
# TODO: You may want to add flags here?
CXXFLAGS += -DNDEBUG
###############################################################################
else
CXXFLAGS += -O0 -g
endif

###############################################################################
# TODO: Write the main make target 'vec_red'.  This target will compile your
# code using the flags you set in CXXFLAGS.  Be sure that the generated
# executable is named 'vec_red' in order to use the BASH script
# measure_speedup.sh later on.
vec_red:
###############################################################################

# You can submit a job on euler using this target.  It depends on the 'vec_red'
# target, that is, it will compile the code if there are changes prior to
# submitting the job.  It will then submit the contents of the script
# measure_speedup.sh
measurement: vec_red
	bsub -W 00:15 \
		-n 24 -R fullnode -R 'rusage[scratch=16]' < measure_speedup.sh

clean:
	rm -f vec_red
