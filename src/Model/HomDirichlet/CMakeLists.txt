cmake_minimum_required(VERSION 2.8.5)

project(HomDirichlet CXX)

# Include configuration/headers/libraries common to all code
INCLUDE(${QUINOA_SOURCE_DIR}/QuinoaCommon.cmake)


include_directories(${QUINOA_SOURCE_DIR}/Model/HomDirichlet
                    ${QUINOA_SOURCE_DIR}/Model/MixModel/Dirichlet
                    ${QUINOA_SOURCE_DIR}/Model/MixModel
                    ${QUINOA_SOURCE_DIR}/Model
                    ${QUINOA_SOURCE_DIR}/Random
                    ${QUINOA_SOURCE_DIR}/Base
)

add_library(HomDirichlet
            HomDirichlet.C
)

INSTALL(TARGETS HomDirichlet
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
