#!/bin/sh

# At this time this script is not able to read $(MAKE) from the test harness;
# we are assuming 'make' will do the job

cd ../
make serialtests
cd tests-core

