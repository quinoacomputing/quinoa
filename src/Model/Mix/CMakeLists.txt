cmake_minimum_required(VERSION 2.8.5)

project(Mix CXX)

include_directories(${QUINOA_SOURCE_DIR}/Model/Mix
                    ${QUINOA_SOURCE_DIR}/Model
                    ${QUINOA_SOURCE_DIR}/Model/Mix/Dirichlet
                    ${QUINOA_SOURCE_DIR}/Statistics
                    ${QUINOA_SOURCE_DIR}/Control
                    ${QUINOA_SOURCE_DIR}/Random
                    ${QUINOA_SOURCE_DIR}/Base
)

add_library(Mix
            Mix.C
)

INSTALL(TARGETS Mix
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
