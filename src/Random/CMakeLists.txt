cmake_minimum_required(VERSION 2.8.5)

project(Random CXX)

include_directories(${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Control
                    ${QUINOA_SOURCE_DIR}/Paradigm
                    ${QUINOA_SOURCE_DIR}/Random
)

add_library(Random
            Random.C
            MKL.C
            MKLRandom.C
            MKLRndStream.C
            MKLRndTable.C
            MKLTest.C
)

INSTALL(TARGETS Random
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
