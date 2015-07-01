@TEST_MAKE_INCLUDE@
AM_LDFLAGS = @AM_LDFLAGS@
AM_CPPFLAGS = @AM_CPPFLAGS@
LDADD = $(top_builddir)/src/libmesquite.la
check_PROGRAMS = main
main_SOURCES = main.cpp
main_DEPENDENCIES = $(LDADD)
TESTS = $(check_PROGRAMS)
MOSTLYCLEANFILES = *.vtk *.gpt *.gnuplot
EXTRA_DIST = CMakeLists.txt
