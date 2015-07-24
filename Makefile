# |
# o---------------------------------------------------------------------o
# |
# | Numdiff makefile
# |
# o---------------------------------------------------------------------o
# |
# | Methodical Accelerator Design
# |
# | Copyright (c) 2011+ CERN, mad@cern.ch
# |
# | For more information, see http://cern.ch/mad/numdiff
# |
# o---------------------------------------------------------------------o
# |
# | $Id$
# |

# For makefile documentation, please read make/README
# For information and bug report, please contact mad@cern.ch
#

###################
# Project settings

PROJECT := numdiff

#################
# Build settings
#

# architecture bit: detect/32/64 (default is detect)
ARCH    := detect

# debugging mode: yes/no (default is no)
DEBUG   := no

# profiling mode: yes/no (default is no)
PROFILE := no

#############################
# Compilers/Linkers settings
# see make/compiler.* for supported compilers
# GNU=yes   sets CC=gcc,     CXX=g++,     FC=gfortran (default)
# Intel=yes sets CC=icc/icl, CXX=icc/icl, FC=ifort    (use icl on Windows)

# C compiler (default is gcc)
CC := gcc

# Linker (default is C compiler)
LD  = $(CC)

# Tester (default is numdiff)
ND  := numdiff

####################
# Makefile includes

makedir := ../../make
include $(makedir)/make.inc

# end of makefile
