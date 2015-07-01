@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/src/libmesquite.la
check_PROGRAMS = par_hex_untangle_shape
par_hex_untangle_shape_SOURCES = par_hex_untangle_shape.cpp
par_hex_untangle_shape_DEPENDENCIES = $(LDADD)
TESTS = $(check_PROGRAMS)
EXTRA_DIST = CMakeLists.txt

MOSTLYCLEANFILES = *.*.vtk *.gpt
