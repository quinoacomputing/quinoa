cmake_minimum_required(VERSION 2.8.5)

project(Geometry CXX)

include_directories(${QUINOA_SOURCE_DIR}/Geometry
                    ${QUINOA_SOURCE_DIR}/Control
                    ${QUINOA_SOURCE_DIR}/Mesh
                    ${QUINOA_SOURCE_DIR}/IO
                    ${QUINOA_SOURCE_DIR}/Base
                    ${TPL_DIR}/include
)

add_library(Geometry
            Geometry.C
            AnalyticGeometry.C
            DiscreteGeometry.C
)

INSTALL(TARGETS Geometry
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
