cmake_minimum_required(VERSION 2.8.5)

project(WallTreatment CXX)

include_directories(${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/LinearAlgebra
)

add_library(WallTreatment
            EllipticRelaxation.C
)

INSTALL(TARGETS WallTreatment
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
