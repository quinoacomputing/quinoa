if ENABLE_IMESH
  IMESH_LIB = libmsqimesh.la
else
  IMESH_LIB = 
endif

if ENABLE_IGEOM
  IGEOM_LIB = libmsqigeom.la
else
  IGEOM_LIB = 
endif

if ENABLE_IREL
  IREL_LIB = libmsqirel.la
else
  IREL_LIB = 
endif

if ENABLE_IMESHP
  IMESHP_LIB = libmsqimeshp.la
else
  IMESHP_LIB = 
endif


lib_LTLIBRARIES = $(IMESH_LIB) $(IGEOM_LIB) $(IREL_LIB) $(IMESHP_LIB)

libmsqimesh_la_includedir = $(includedir)
libmsqimesh_la_SOURCES = MsqIMesh.cpp 
libmsqimesh_la_include_HEADERS = MsqIMesh.hpp 
libmsqimesh_la_LIBADD = $(top_builddir)/src/libmesquite.la $(IMESH_LIBS)

libmsqigeom_la_includedir = $(includedir)
libmsqigeom_la_SOURCES = MsqIGeom.cpp 
libmsqigeom_la_include_HEADERS = MsqIGeom.hpp 
libmsqigeom_la_LIBADD = $(top_builddir)/src/libmesquite.la $(IGEOM_LIBS)

libmsqirel_la_includedir = $(includedir)
libmsqirel_la_SOURCES = MsqIRel.cpp 
libmsqirel_la_include_HEADERS = MsqIRel.hpp 
libmsqirel_la_LIBADD = libmsqigeom.la libmsqimesh.la $(IREL_LIBS)

libmsqimeshp_la_includedir = $(includedir)
libmsqimeshp_la_SOURCES = MsqIMeshP.cpp 
libmsqimeshp_la_include_HEADERS = MsqIMeshP.hpp 
libmsqimeshp_la_LIBADD = libmsqimesh.la $(IMESHP_LIBS)

AM_CPPFLAGS = $(CPPFLAGS) \
              -I$(top_builddir)/src/include \
              -I$(top_srcdir)/src/include \
              -I$(top_srcdir)/src/Mesh \
              -I$(top_srcdir)/src/Misc \
              -I$(top_srcdir)/src/Control \
              -I$(top_srcdir)/src/MappingFunction \
              $(IMESH_INCL) \
              $(IGEOM_INCL) \
              $(IREL_INCL) \
              $(IMESHP_INCL) 

EXTRA_DIST = CMakeLists.txt



