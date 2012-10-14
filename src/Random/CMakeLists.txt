cmake_minimum_required(VERSION 2.8.5)

project(Random CXX)

# Include configuration/headers/libraries common to all code
INCLUDE(${QUINOA_SOURCE_DIR}/QuinoaCommon.cmake)


include_directories(${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Random
)

add_library(Random
            Random.C
            RandomException.C
            MKLException.C
            MKLRandom.C
)

INSTALL(TARGETS Random
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
