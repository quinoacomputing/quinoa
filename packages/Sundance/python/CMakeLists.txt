


# Simple wrapper interface generator (SWIG) is required
FIND_PACKAGE(SWIG 1.3.29 REQUIRED)
# ${SWIG_USE_FILE} is hard-coded to be the installed UseSWIG.cmake
# file.  We want to include the local, modified version.
INCLUDE(UseSWIG)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(example)


