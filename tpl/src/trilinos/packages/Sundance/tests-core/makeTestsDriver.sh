#!/bin/sh

# At this time this script is not able to read $(MAKE) from the test harness;
# we are assuming 'make' will do the job
# The current setup also does not recognize proc maxes and mpigo from the 
# test harness.

cd ../
make tests
cd tests-core

#make tests
#cd ../tests-solvers
#make tests
#cd ../tests-std-framework
#make tests
#cd ../tests-std-mesh
#make tests
#cd ../tests-utils
#make tests
#cd ../tests-core
