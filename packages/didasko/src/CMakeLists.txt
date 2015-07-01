
#This file is only here to set up the config file in the right spot.
#potentially we could put the latex building logic in here as well at
#some point.

INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

IF(HAVE_DIDASKO_EPETRAEXT AND ${PACKAGE_NAME}_ENABLE_Zoltan)
  SET(HAVE_ZOLTAN ON)
ENDIF()

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
