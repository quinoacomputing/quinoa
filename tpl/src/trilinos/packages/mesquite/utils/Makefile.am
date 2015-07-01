@UTIL_MAKE_INCLUDE@
LDADD = libmsqutil.la $(top_builddir)/src/libmesquite.la
bin_PROGRAMS = vtkrandom vtkxform msqshape msqquality
man1_MANS = vtkrandom.1 vtkxform.1 msqshape.1 msqquality.1
noinst_LTLIBRARIES = libmsqutil.la

libmsqutil_la_SOURCES = CLArgs.cpp CLArgs.hpp \
                        CLArgFlag.cpp CLArgFlag.hpp \
                        ManPage.cpp ManPage.hpp

vtkrandom_SOURCES = randomize.cpp domain.hpp domain.cpp
msqshape_SOURCES = shape-improve.cpp domain.hpp domain.cpp
vtkxform_SOURCES = transform.cpp
msqquality_SOURCES = quality.cpp domain.hpp domain.cpp

vtkrandom.1: vtkrandom$(EXEEXT)
	./$< -M > $@ || touch $@

msqshape.1: msqshape$(EXEEXT)
	./$< -M > $@ || touch $@

vtkxform.1: vtkxform$(EXEEXT)
	./$< -M > $@ || touch $@

msqquality.1: msqquality$(EXEEXT)
	./$< -M > $@ || touch $@

MOSTLYCLEANFILES = $(man1_MANS)

EXTRA_DIST = CMakeLists.txt
