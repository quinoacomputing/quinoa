#!/bin/bash
# Return the number of CPUs
echo `cat /proc/cpuinfo | grep MHz | wc -l`
