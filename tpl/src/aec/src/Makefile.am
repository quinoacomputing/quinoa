lib_LTLIBRARIES = libaec.la libsz.la
libaec_la_SOURCES = encode.c encode_accessors.c decode.c \
encode.h encode_accessors.h decode.h
libaec_la_LDFLAGS = -version-info 0:3:0

libsz_la_SOURCES = sz_compat.c
libsz_la_LIBADD = libaec.la
libsz_la_LDFLAGS = -version-info 2:1:0

include_HEADERS = libaec.h szlib.h

bin_PROGRAMS = aec
noinst_PROGRAMS = utime
utime_SOURCES = utime.c
aec_LDADD = libaec.la
aec_SOURCES = aec.c

EXTRA_DIST = CMakeLists.txt
