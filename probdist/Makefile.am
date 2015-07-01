include $(top_srcdir)/include/Makefile.def

lib_LTLIBRARIES = libprobdist.la
libprobdist_la_SOURCES = $(PROBDISTSOURCES)
libprobdist_la_LIBADD = $(top_builddir)/mylib/libmylib.la
libprobdist_la_LDFLAGS = -no-undefined -version-info 0:1:0

EXTRA_DIST = guideprobdist.tex titre.tex copyright.tex guideprobdist.bbl

doc_DATA = guideprobdist.pdf
