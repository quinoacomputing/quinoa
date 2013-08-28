cmake_minimum_required(VERSION 2.8.5)

project(GeneralizedDirichlet CXX)

include_directories(${QUINOA_SOURCE_DIR}/Model/Mix/GeneralizedDirichlet
                    ${QUINOA_SOURCE_DIR}/Model/Mix
                    ${QUINOA_SOURCE_DIR}/Model
                    ${QUINOA_SOURCE_DIR}/Random
                    ${QUINOA_SOURCE_DIR}/Statistics
                    ${QUINOA_SOURCE_DIR}/Control
                    ${QUINOA_SOURCE_DIR}/Base
)

add_library(GeneralizedDirichlet
            GeneralizedDirichlet.C
)

INSTALL(TARGETS GeneralizedDirichlet
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
