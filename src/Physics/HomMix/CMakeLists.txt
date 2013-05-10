cmake_minimum_required(VERSION 2.8.5)

project(HomMix CXX)

# Include configuration/headers/libraries common to all code
INCLUDE(${QUINOA_SOURCE_DIR}/QuinoaCommon.cmake)


include_directories(${QUINOA_SOURCE_DIR}/Physics/HomMix
                    ${QUINOA_SOURCE_DIR}/Model
                    ${QUINOA_SOURCE_DIR}/Model/Mix
                    ${QUINOA_SOURCE_DIR}/Model/Hydro
                    ${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Physics
                    ${QUINOA_SOURCE_DIR}/Statistics
                    ${QUINOA_SOURCE_DIR}/IO
                    ${QUINOA_SOURCE_DIR}/Control
)

add_library(HomMix
            HomMix.C
)

INSTALL(TARGETS HomMix
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
