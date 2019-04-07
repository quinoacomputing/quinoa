################################################################################
#
# \file      FindMCSS.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find m.css
#
################################################################################

# m.css: http://mcss.mosra.cz
#
#  MCSS_FOUND - System has numdiff
#  MCSS_DOX2HTML5 - The dox2html5 python script
#
#  Usage:
#
#  set(MCSS_ROOT "/path/to/custom/mcss") # prefer over system
#  find_package(MCSS)

if(MCSS_DOX2HTML5 AND PYTHONINTERP_FOUND AND PYGMENTS_FOUND AND JINJA2_FOUND AND LATEX_FOUND)
  # Already in cache, be silent
  set (MCSS_FIND_QUIETLY TRUE)
endif()

# Find Python 3 and prerequisites for m.css
find_package(PythonInterp 3.6)

if (PYTHONINTERP_FOUND)

  execute_process( COMMAND ${PYTHON_EXECUTABLE} -c "import pygments"
                   ERROR_VARIABLE PYGMENTS_STDERR)
  execute_process( COMMAND ${PYTHON_EXECUTABLE} -c "import jinja2"
                   ERROR_VARIABLE JINJA2_STDERR)

  if (PYGMENTS_STDERR)
    set(PYGMENTS_FOUND "false")
    message(STATUS "Could NOT find Python module Pygments")
    #message("${PYGMENTS_STDERR}")
  else()
    set(PYGMENTS_FOUND "true")
    message(STATUS "Found Python module Pygments")
  endif()

  if (JINJA2_STDERR)
    set(JINJA2_FOUND "false")
    message(STATUS "Could NOT find Python module Jinja2")
    #message("${JINJA2_STDERR}")
  else()
    set(JINJA2_FOUND "true")
    message(STATUS "Found Python module Jinja2")
  endif()

endif()

find_package(LATEX)

FIND_PROGRAM(MCSS_DOX2HTML5 NAMES dox2html5.py
                            PATHS ${MCSS_ROOT} $ENV{MCSS_ROOT}
                            PATH_SUFFIXES m.css/doxygen)

# Handle the QUIETLY and REQUIRED arguments and set MCSS_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MCSS DEFAULT_MSG MCSS_DOX2HTML5 PYTHONINTERP_FOUND PYGMENTS_FOUND JINJA2_FOUND LATEX_FOUND)

MARK_AS_ADVANCED(MCSS_DOX2HTML5 PYTHONINTERP_FOUND PYGMENTS_FOUND JINJA2_FOUND LATEX_FOUND)
