@TEST_MAKE_INCLUDE@
AM_LDFLAGS = @AM_LDFLAGS@
AM_CPPFLAGS = @AM_CPPFLAGS@
LDADD = $(top_builddir)/src/libmesquite.la
AM_CPPFLAGS += '-DSRCDIR="$(srcdir)/"'
check_PROGRAMS = main
main_SOURCES = main.cpp
main_DEPENDENCIES = $(LDADD)
TESTS = $(check_PROGRAMS)
EXTRA_DIST = CMakeLists.txt reference.vtk figures.tex
MOSTLYCLEANFILES = *.gpt Exp*.vtk
