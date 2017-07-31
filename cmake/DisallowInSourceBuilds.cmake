################################################################################
#
# \file      cmake/DisallowInSourceBuilds.cmake
# \copyright 2016-2017, Los Alamos National Security, LLC.
# \brief     Cmake code to disallow in-source builds
#
################################################################################

if(__disallow_in_source_builds)
	return()
endif()
set(disallow_in_source_builds YES)


function(disallow_in_source_builds)

  # Disallow in-source builds
  if("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
    message(FATAL_ERROR "ERROR! In-source builds are not allowed!\n"
    "Please create a separate directory and configure there; for example:\n"
    "  $ mkdir build\n"
    "  $ cd build\n"
    "  $ cmake [OPTIONS] ..."
    "\n"
    "NOTE: You must first delete the newly created CMakeCache.txt file and"
    "CMakeFiles directory under the source directory or you will not be able to"
    "configure correctly! For example:\n"
    "  $ rm -r CMakeCache.txt CMakeFiles")
  endif()

endfunction()
