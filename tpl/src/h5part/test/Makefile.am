# test level Makefile.am

OBJEXT=o

INCLUDES = -I../src @INCLUDES@

noinst_PROGRAMS = @BUILD_TESTS@

TESTS_ENVIRONMENT = env LD_LIBRARY_PATH=@HDF5ROOT@/lib:$(LD_LIBRARY_PATH)
TESTS = @BUILD_TESTS@

EXTRA_PROGRAMS = test testf

test_SOURCES =	test.c \
		testframe.c \
		write.c \
		read.c \
		testframe.h \
		params.h
test_LDADD   =  ../src/libH5Part.la

testf_SOURCES = testf.F90
testf_LDADD   =  ../src/libH5PartF.la

clean-local:
	$(RM) -f *.h5

