@TEST_MAKE_INCLUDE@
AM_LDFLAGS = @AM_LDFLAGS@
AM_CPPFLAGS = @AM_CPPFLAGS@
LDADD = $(top_builddir)/src/libmesquite.la
AM_CPPFLAGS += '-DSRCDIR="$(srcdir)/"'
check_PROGRAMS = main 
main_SOURCES = main.cpp
main_DEPENDENCIES = $(LDADD)
TESTS = $(check_PROGRAMS)
EXTRA_DIST = CMakeLists.txt initial.vtk 
MOSTLYCLEANFILES = 1-1.results.vtk \
                   1-2.results.vtk \
                   1-4.results.vtk \
                   1-8.results.vtk \
                   2-NE.results.vtk \
                   2-NW.results.vtk \
                   2-SE.results.vtk \
                   2-SW.results.vtk \
                   3-Left.results.vtk \
                   3-Mid.results.vtk \
                   3-Right.results.vtk

