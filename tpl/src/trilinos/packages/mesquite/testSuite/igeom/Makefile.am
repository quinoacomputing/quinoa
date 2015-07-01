@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/itaps/libmsqigeom.la \
        $(top_builddir)/src/libmesquite.la \
        $(IGEOM_LIBS)
check_PROGRAMS = main
main_SOURCES = main.cpp
main_DEPENDENCIES = $(top_builddir)/src/libmesquite.la \
                    $(top_builddir)/itaps/libmsqigeom.la
TESTS = main
EXTRA_DIST = CMakeLists.txt 

AM_CPPFLAGS += $(IGEOM_INCL)

if ENABLE_IGEOM

else
main$(EXEEXT):
	echo "ITAPS support disabled.  Reconfigure to build this test."
endif

MOSTLYCLEANFILES =	
