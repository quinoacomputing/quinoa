cmake_minimum_required(VERSION 2.8.5)

project(Random CXX)

# Include configuration/headers/libraries common to all code
INCLUDE(${QUINOA_SOURCE_DIR}/QuinoaCommon.cmake)


include_directories(${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Random
                    ${QUINOA_SOURCE_DIR}/Paradigm
)

add_library(Random
            Random.C
            RandomException.C
            MKL.C
            MKLRandom.C
            MKLException.C
            VSLException.C
            MKLRndStream.C
            MKLRndTable.C
)

INSTALL(TARGETS Random
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
