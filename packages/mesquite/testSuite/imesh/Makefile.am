@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/itaps/libmsqimesh.la \
        $(top_builddir)/src/libmesquite.la \
        $(IMESH_LIBS)
check_PROGRAMS = main
main_SOURCES = main.cpp
main_DEPENDENCIES = $(top_builddir)/src/libmesquite.la \
                    $(top_builddir)/itaps/libmsqimesh.la
TESTS = main
EXTRA_DIST = CMakeLists.txt

AM_CPPFLAGS += $(IMESH_INCL)

if ENABLE_IMESH

else
main$(EXEEXT):
	echo "ITAPS support disabled.  Reconfigure to build this test."
endif

MOSTLYCLEANFILES =	feasible-newton-result.vtk \
			original.vtk \
			smart-laplacian-result.vtk
