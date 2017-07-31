################################################################################
#
# \file      cmake/ConfigureDataLayout.cmake
# \copyright 2016-2017, Los Alamos National Security, LLC.
# \brief     Configure data layouts
#
################################################################################

# Configure data layout for particle data

# Available options
set(PARTICLE_DATA_LAYOUT_VALUES "particle" "equation")
# Initialize all to off
set(PARTICLE_DATA_LAYOUT_AS_PARTICLE_MAJOR off)  # 0
set(PARTICLE_DATA_LAYOUT_AS_EQUATION_MAJOR off)  # 1
# Set default and select from list
set(PARTICLE_DATA_LAYOUT "particle" CACHE STRING "Particle data layout. Default: (particle-major). Available options: ${PARTICLE_DATA_LAYOUT_VALUES}(-major).")
SET_PROPERTY (CACHE PARTICLE_DATA_LAYOUT PROPERTY STRINGS ${PARTICLE_DATA_LAYOUT_VALUES})
STRING (TOLOWER ${PARTICLE_DATA_LAYOUT} PARTICLE_DATA_LAYOUT)
LIST (FIND PARTICLE_DATA_LAYOUT_VALUES ${PARTICLE_DATA_LAYOUT} PARTICLE_DATA_LAYOUT_INDEX)
# Evaluate selected option and put in a define for it
IF (${PARTICLE_DATA_LAYOUT_INDEX} EQUAL 0)
  set(PARTICLE_DATA_LAYOUT_AS_PARTICLE_MAJOR on)
ELSEIF (${PARTICLE_DATA_LAYOUT_INDEX} EQUAL 1)
  set(PARTICLE_DATA_LAYOUT_AS_EQUATION_MAJOR on)
ELSEIF (${PARTICLE_DATA_LAYOUT_INDEX} EQUAL -1)
  MESSAGE(FATAL_ERROR "Particle data layout '${PARTICLE_DATA_LAYOUT}' not supported, valid entries are ${PARTICLE_DATA_LAYOUT_VALUES}(-major).")
ENDIF()
message(STATUS "Particle data layout: " ${PARTICLE_DATA_LAYOUT} "(-major)")

# Configure data layout for mesh field data

# Available options
set(FIELD_DATA_LAYOUT_VALUES "field" "equation")
# Initialize all to off
set(FIELD_DATA_LAYOUT_AS_FIELD_MAJOR off)  # 0
set(FIELD_DATA_LAYOUT_AS_EQUATION_MAJOR off)  # 1
# Set default and select from list
set(FIELD_DATA_LAYOUT "field" CACHE STRING "Mesh field data layout. Default: (field-major). Available options: ${FIELD_DATA_LAYOUT_VALUES}(-major).")
SET_PROPERTY (CACHE FIELD_DATA_LAYOUT PROPERTY STRINGS ${FIELD_DATA_LAYOUT_VALUES})
STRING (TOLOWER ${FIELD_DATA_LAYOUT} FIELD_DATA_LAYOUT)
LIST (FIND FIELD_DATA_LAYOUT_VALUES ${FIELD_DATA_LAYOUT} FIELD_DATA_LAYOUT_INDEX)
# Evaluate selected option and put in a define for it
IF (${FIELD_DATA_LAYOUT_INDEX} EQUAL 0)
  set(FIELD_DATA_LAYOUT_AS_FIELD_MAJOR on)
ELSEIF (${FIELD_DATA_LAYOUT_INDEX} EQUAL 1)
  set(FIELD_DATA_LAYOUT_AS_EQUATION_MAJOR on)
ELSEIF (${FIELD_DATA_LAYOUT_INDEX} EQUAL -1)
  MESSAGE(FATAL_ERROR "Mesh field data layout '${FIELD_DATA_LAYOUT}' not supported, valid entries are ${FIELD_DATA_LAYOUT_VALUES}(-major).")
ENDIF()
message(STATUS "Mesh field data layout: " ${FIELD_DATA_LAYOUT} "(-major)")
