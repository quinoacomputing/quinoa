# tools level Makefile.am

OBJEXT=o

INCLUDES = -I../src @INCLUDES@

# What to build... make install will place these files in the $(prefix)/bin directory.
bin_PROGRAMS = @BUILD_TOOLS@

# Listing of all programs that maybe built. (Has to know statically...)
EXTRA_PROGRAMS =	h5pAttrib \
			h5pToGNUplot \
			homdynToH5p \
			H5PartBench \
			H5BlockBench

# Listing of sources
h5pAttrib_SOURCES =	h5pAttrib.cc
h5pAttrib_LDADD   =  ../src/libH5Part.la

h5pToGNUplot_SOURCES =	h5pToGNUplot.cc
h5pToGNUplot_LDADD   =  ../src/libH5Part.la

homdynToH5p_SOURCES =	homdynToH5p.cc
homdynToH5p_LDADD   =  ../src/libH5Part.la

H5PartBench_SOURCES =	H5PartBench.c
H5PartBench_LDADD   =  ../src/libH5Part.la

H5BlockBench_SOURCES =	H5BlockBench.c
H5BlockBench_LDADD   =  ../src/libH5Part.la

