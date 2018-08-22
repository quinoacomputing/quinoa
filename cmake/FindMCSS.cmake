################################################################################
#
# \file      cmake/FindMCSS.cmake
# \copyright 2016-2018, Los Alamos National Security, LLC.
# \brief     Find m.css
#
################################################################################

# Find m.css (http://mcss.mosra.cz)
#
#  MCSS_FOUND - System has numdiff
#  MCSS_DOX2HTML5 - The dox2html5 python script
#
#  Usage:
#
#  set(MCSS_ROOT "/path/to/custom/mcss") # prefer over system
#  find_package(MCSS)

if(MCSS_DOX2HTML5)
  # Already in cache, be silent
  set (MCSS_FIND_QUIETLY TRUE)
endif()

find_package(PythonInterp 3.0)

FIND_PROGRAM(MCSS_DOX2HTML5 NAMES dox2html5.py
                            PATHS ${MCSS_ROOT} $ENV{MCSS_ROOT}
                            PATH_SUFFIXES m.css/doxygen)

# Handle the QUIETLY and REQUIRED arguments and set MCSS_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MCSS DEFAULT_MSG MCSS_DOX2HTML5)

get_filename_component(${MCSS_DOX2HTML5} MCSS_DOX2HTML5 ABSOLUTE)
MARK_AS_ADVANCED(MCSS_DOX2HTML5)
