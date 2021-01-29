################################################################################
#
# \file      libstdcxx.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2020 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find libc++ and offer to switch between libc++ and libstdc++
#
################################################################################

if(__libstdcxx)
	return()
endif()
set(__libstdcxx YES)

#### Attempt to find the libc++ library. Do not offer for gnu.
# More infor: libstdc++: http://gcc.gnu.org/libstdc++,
# libc++: http://libcxx.llvm.org.
if (NOT NO_SYSTEM_LIBCXX AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  find_package(LIBCXX)
endif()

#### If libc++ found, offer switch between libstdc++ and libc++.
if (LIBCXX_FOUND)

  # Offer switch between libstdc++ and libc++ for the C++ standard library, only
  # if libc++ has been found. We assume libstdc++ exists.
  function(setlibstdcxx LIBCXX_DEFAULT)

    # Available options
    set(LIBCXX_VALUES "libstdc++" "libc++")

    # Initialize all to off
    set(LIBCXX_STDCPP off)  # 0
    set(LIBCXX_CPP off)     # 1
    # Set default and select from list
    set(LIBCXX ${LIBCXX_DEFAULT} CACHE STRING "Select standard C++ library. Available options: ${LIBCXX_VALUES}.")
    SET_PROPERTY (CACHE LIBCXX PROPERTY STRINGS ${LIBCXX_VALUES})
    STRING (TOLOWER ${LIBCXX} LIBCXX)
    LIST (FIND LIBCXX_VALUES ${LIBCXX} LIBCXX_INDEX)
    # Evaluate selected option and put in a define for it
    IF (${LIBCXX_INDEX} EQUAL 0)
      set(LIBCXX_STDCPP on PARENT_SCOPE)
    ELSEIF (${LIBCXX_INDEX} EQUAL 1)
      set(LIBCXX_CPP on PARENT_SCOPE)
    ELSEIF (${LIBCXX_INDEX} EQUAL -1)
      MESSAGE(FATAL_ERROR "Standard C++ library '${LIBCXX}' not supported, valid entries are ${LIBCXX_VALUES}.")
    ENDIF()
    message(STATUS "Standard C++ library (LIBCXX): " ${LIBCXX})

  endfunction()

  # Evaluate offer for libstdc++ vs libc++ and set compiler flags
  # Arguments: FLAGS2ADD2 - The cmake variable to which to add libc++ args to
  #            DEFAULT    - The default to use if user did not select anything
  macro(set_libstdcpp_vs_libcpp FLAGS2ADD2 DEFAULT)
    if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      if (NOT STDLIBCPP)
        setlibstdcxx(${DEFAULT})
      else()
        setlibstdcxx(${STDLIBCPP})
      endif()
      if(LIBCXX_CPP)  # set compiler flags for use of libc++
        set(${FLAGS2ADD2} "${${FLAGS2ADD2}} -stdlib=libc++")
      else()  # set compiler flags for use of libstdc++
        set(${FLAGS2ADD2} "${${FLAGS2ADD2}} -stdlib=libstdc++")
      endif()
    endif()
  endmacro()

endif()
