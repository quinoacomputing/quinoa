@TEST_MAKE_INCLUDE@
AM_LDFLAGS = @AM_LDFLAGS@
AM_CPPFLAGS = @AM_CPPFLAGS@
AM_CPPFLAGS += '-DSRCDIR="$(srcdir)/"'
LDADD = $(top_builddir)/src/libmesquite.la
check_PROGRAMS = main voshell
TESTS = $(check_PROGRAMS)

main_SOURCES = main.cpp
voshell_SOURCES = voshell.cpp

EXTRA_DIST = CMakeLists.txt \
             expected_linear_output.vtk \
             expected_quadratic_output.vtk \
             hour-quad8.vtk \
	     linear_input.vtk \
	     quadratic_input.vtk \
	     cutCube.vtk \
	     hex3Dpart.vtk \
	     homogeneousPart.vtk \
	     quadPart.vtk \
	     sphereCube.vtk \
	     sphereCylinder_1194_inv.vtk \
	     triPart.vtk

MOSTLYCLEANFILES = smooth_ho.vtk
