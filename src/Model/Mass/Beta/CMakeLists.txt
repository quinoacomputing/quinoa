cmake_minimum_required(VERSION 2.8.5)

project(Beta CXX)

include_directories(${QUINOA_SOURCE_DIR}/Model/Mass/Beta
                    ${QUINOA_SOURCE_DIR}/Model/Mass
                    ${QUINOA_SOURCE_DIR}/Model
                    ${QUINOA_SOURCE_DIR}/Random
                    ${QUINOA_SOURCE_DIR}/Statistics
                    ${QUINOA_SOURCE_DIR}/Control
                    ${QUINOA_SOURCE_DIR}/Base
)

add_library(Beta
            Beta.C
)

INSTALL(TARGETS Beta
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
