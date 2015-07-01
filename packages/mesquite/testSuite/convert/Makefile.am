@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/src/libmesquite.la
check_PROGRAMS = exo_to_vtk vtk_to_exo
exo_to_vtk_SOURCES = exo_to_vtk.cpp
vtk_to_exo_SOURCES = vtk_to_exo.cpp
exo_to_vtk_DEPENDENCIES = $(LDADD)
vtk_to_exo_DEPENDENCIES = $(LDADD)
EXTRA_DIST = CMakeLists.txt
