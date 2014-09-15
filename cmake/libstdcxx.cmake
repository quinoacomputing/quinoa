if(__libstdcxx)
	return()
endif()
set(__libstdcxx YES)


# Set C++ standard library; this should be the same for the third-party
# libraries as well as all executables

function(setlibstdcxx LIBCXX_DEFAULT)

  # Available options, for more info, see:
  # libstdc++: http://gcc.gnu.org/libstdc++
  # libc++: http://libcxx.llvm.org
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
