cmake_minimum_required(VERSION 2.8.5)

project(HomogeneousGeneralizedDirichlet CXX)

# Include configuration/headers/libraries common to all code
INCLUDE(${QUINOA_SOURCE_DIR}/QuinoaCommon.cmake)


include_directories(${QUINOA_SOURCE_DIR}/Physics/HomogeneousGeneralizedDirichlet
                    ${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Random
                    ${QUINOA_SOURCE_DIR}/Physics
)

add_library(HomogeneousGeneralizedDirichlet
            HomogeneousGeneralizedDirichlet.C
)

INSTALL(TARGETS HomogeneousGeneralizedDirichlet
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
