
INCLUDE(TribitsPackageMacros)

TRIBITS_PACKAGE(ForTrilinos Fortran)

IF("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qfree=f90 -qxlf2003=polymorphic -qxlf2003=autorealloc -WF,-C!")
#    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qxlf2003=polymorphic -qxlf2003=autorealloc -qxlf90=autodealloc -WF,-C!")
    SET(ForTrilinos_UT_CPP_FLAGS -E -P -qnoppline)
    SET(ForTrilinos_UT_CPP ${CMAKE_C_COMPILER})
ELSE()
    SET(ForTrilinos_UT_CPP_FLAGS "")
    SET(ForTrilinos_UT_CPP cpp)
ENDIF()

# Previously, unit testing required support for long lines
# Enabled by default for: "NAG", "XL", Intel
# Limited support: "PGI" (limited to ~250)
# Flag needed for: "GNU"
IF("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
ENDIF()
#IF("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL|NAG|Intel|GNU")
#    SET(FORTRILINOS_COMPILER_SUPPORTS_LONG_LINES 1)
#ENDIF()
# NOT: Now, unit testing only requires Perl be installed and in the path

# A) Package-specific configuration options
#

# Normalize flags for preprocessing by Fortran compilers
FOREACH(var
    HAVE_FORTRILINOS_EXAMPLES
    HAVE_FORTRILINOS_TESTS
    HAVE_EXAMPLES
    HAVE_EXPORT_MAKEFILES
    HAVE_GNUMAKE
    ForTrilinos_ASSERTIONS
    HAVE_MPI
    HAVE_TESTS
    HAVE_TEUCHOS_COMPLEX
    ForTrilinos_ENABLE_EXTENDED  # DEPRECATED Flag remove in Trilinos 11.0, change to ForTrilinos_ENABLE_OBJECT_ORIENTED
    ForTrilinos_ENABLE_OBJECT_ORIENTED
    ForTrilinos_DISABLE_DEFERRED_LENGTH_CHARACTERS
    ForTrilinos_DISABLE_FINAL_SUBROUTINES
    HAVE_FORTRILINOS_AMESOS
    ForTrilinos_ENABLE_Amesos
    HAVE_FORTRILINOS_AZTECOO
    ForTrilinos_ENABLE_AztecOO
    HAVE_FORTRILINOS_IFPACK
    ForTrilinos_ENABLE_Ifpack
    HAVE_FORTRILINOS_GALERI
    ForTrilinos_ENABLE_Galeri
    HAVE_FORTRILINOS_PLIRIS
    ForTrilinos_ENABLE_Pliris
    )
  IF(${var})
    SET(ForTrilinos_config_${var} 1)
  ELSE()
    SET(ForTrilinos_config_${var} 0)
  ENDIF()
ENDFOREACH()

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)




ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

TRIBITS_PACKAGE_POSTPROCESS()
