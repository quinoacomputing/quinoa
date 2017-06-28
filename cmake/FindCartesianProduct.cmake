################################################################################
#
# \file      cmake/FindCartesianProduct.cmake
# \copyright 2016, Los Alamos National Security, LLC.
# \brief     Find the Cartesian Product header
#
################################################################################

# Find Boost MPL CartesianProduct headers
#
#  CARTESIAN_PRODUCT_FOUND        - True if Cartesian Product is found
#  CARTESIAN_PRODUCT_INCLUDE_DIRS - Cartesian Product include files paths
#
#  Set CARTESIAN_PRODUCT_ROOT before calling find_package to a path to
#  add an additional search path, e.g.,
#
#  Usage:
#
#  set(CARTESIAN_PRODUCT_ROOT "/path/to/custom/boost/mpl") # prefer over system
#  find_package(CartesianProduct)
#  if(CARTESIAN_PRODUCT_FOUND)
#    include_directories(${CARTESIAN_PRODUCT_INCLUDE_DIRS})
#  endif()

# If already in cache, be silent
if (CARTESIAN_PRODUCT_INCLUDE_DIRS)
  set (BoostMPLCartesianProduct_FIND_QUIETLY TRUE)
endif()

# Look for the header file
FIND_PATH(CARTESIAN_PRODUCT_INCLUDE_DIR
          NAMES boost/mpl/cartesian_product.hpp
          HINTS ${CARTESIAN_PRODUCT_ROOT}/include
                $ENV{CARTESIAN_PRODUCT_ROOT})

set(CARTESIAN_PRODUCT_INCLUDE_DIRS ${CARTESIAN_PRODUCT_INCLUDE_DIR})

# Handle the QUIETLY and REQUIRED arguments and set CARTESIAN_PRODUCT_FOUND to
# TRUE if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CartesianProduct DEFAULT_MSG
                                  CARTESIAN_PRODUCT_INCLUDE_DIRS)

MARK_AS_ADVANCED(CARTESIAN_PRODUCT_INCLUDE_DIRS)
