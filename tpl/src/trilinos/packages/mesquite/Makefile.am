AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
LIBTOOL = @LIBTOOL@ @MSQ_LIBTOOL_FLAGS@
AM_LDFLAGS = @AM_LDFLAGS@
AM_CPPFLAGS = @AM_CPPFLAGS@

SUBDIRS = . src itaps doc testSuite utils meshFiles


# Need to be explicit about these for old versions of automake
OLD_AM_EXTRA_CONFIG = m4/cxxflags.m4 \
                      m4/cxx.m4 \
                      m4/util.m4 \
		      m4/itaps.m4
 
EXTRA_DIST = lgpl.txt \
             mkconfigure \
             $(OLD_AM_EXTRA_CONFIG) \
             $(INCLUDED_MAKEFILES) \
             $(srcdir)/cmake/*.cmake \
	     CMakeLists.txt \
             ReleaseNotes.txt

