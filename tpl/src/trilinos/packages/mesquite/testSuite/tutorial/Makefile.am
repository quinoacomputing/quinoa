@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/src/libmesquite.la
check_PROGRAMS = tutorial wrapper lowlevel
TESTS = run
EXTRA_DIST = CMakeLists.txt run

tutorial_SOURCES = tutorial.cpp
tutorial_DEPENDENCIES = $(LDADD)

wrapper_SOURCES = wrapper.cpp
wrapper_DEPENDENCIES = $(LDADD)

lowlevel_SOURCES = lowlevel.cpp
lowlevel_DEPENDENCIES = $(LDADD)

MOSTLYCLEANFILES = original_mesh.vtk  smoothed_mesh.vtk

exampledir = ${docdir}/examples
example_DATA = tutorial.cpp tutorial.make
