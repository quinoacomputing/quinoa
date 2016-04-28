ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src tests
EXTRA_DIST = doc/license.txt CMakeLists.txt cmake/config.h.in \
cmake/macros.cmake README.SZIP

bench: all
	cd tests && $(MAKE) $(AM_MAKEFLAGS) bench
benc: all
	cd tests && $(MAKE) $(AM_MAKEFLAGS) benc
bdec: all
	cd tests && $(MAKE) $(AM_MAKEFLAGS) bdec
.PHONY: bench benc bdec
