# src level Makefile.am

OBJEXT=o

EXTRA_HEADERS = 	H5PartPrivate.h \
			H5BlockPrivate.h \
			H5MultiBlockPrivate.h \
			H5PartAttrib.h

# Extra files that I wish to include in the dist tar ball.
EXTRA_DIST =		TestUnderscoreC.c \
			TestUnderscore.f \
			generate-attr.py \
			generate-h5bl-readwrite.py \
			generate-h5multi-readwrite.py \
			H5Part.f90 \
			H5PartAttrib.f90 \
			H5Block.f90 \
			H5BlockReadWrite.f90 \
			$(EXTRA_HEADERS)

# Files that I don't want to include in the dist tar ball
nodist_include_HEADERS = H5PartF.h @UNDERSCORE_H@

# What to build... Will be determined by configure script.
lib_LTLIBRARIES = libH5Part.la @LIB_FORTRAN@

# Listing of all possible targets that I may build.
EXTRA_LTLIBRARIES = libH5PartF.la

# Header files that I wish to install in $(prefix)/include
include_HEADERS = 	H5Part.h \
			H5PartTypes.h \
			H5PartErrors.h \
			H5PartAttrib.h \
			H5Block.h \
			H5BlockTypes.h \
			H5BlockErrors.h \
			H5BlockReadWrite.h \
			H5MultiBlock.h \
			H5MultiBlockTypes.h \
			H5MultiBlockErrors.h \
			H5MultiBlockReadWrite.h \
			H5PartF.h \
			@UNDERSCORE_H@

# Listing of sources
libH5Part_la_SOURCES =	H5Part.c \
			H5PartAttrib.c \
			H5Block.c \
			H5BlockReadWrite.c \
			H5MultiBlock.c \
			H5MultiBlockReadWrite.c

libH5Part_la_LDFLAGS = -version-info 0:0:0

libH5PartF_la_SOURCES =	H5Part.c \
			H5PartAttrib.c \
			H5Block.c \
			H5BlockReadWrite.c \
			H5MultiBlock.c \
			H5MultiBlockReadWrite.c \
			H5PartF.c \
			H5PartAttribF.c \
			H5BlockF.c \
			H5BlockReadWriteF.c

libH5PartF_la_LDFLAGS = -version-info 0:0:0 -rpath '$(libdir)'

H5PartF.h: H5Part.f90 H5PartAttrib.f90 H5Block.f90 H5BlockReadWrite.f90
	awk '/INTEGER\*8 FUNCTION/{print "\t" $$1 " " $$3}' $^ >$@

