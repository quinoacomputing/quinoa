include $(top_srcdir)/include/Makefile.def

include_HEADERS = TestU01.h

lib_LTLIBRARIES = libtestu01.la
libtestu01_la_SOURCES = $(TESTU01SOURCES)
libtestu01_la_LIBADD = $(top_builddir)/probdist/libprobdist.la $(top_builddir)/mylib/libmylib.la
libtestu01_la_LDFLAGS = -no-undefined -version-info 0:1:0

EXTRA_DIST = guidetestu01.tex intro.tex title.tex scatfig.tex bintro.tex copyright.tex fintro.tex sintro.tex uintro.tex guidetestu01.bbl

doc_DATA = guideshorttestu01.pdf guidelongtestu01.pdf
