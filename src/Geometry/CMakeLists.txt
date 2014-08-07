cmake_minimum_required(VERSION 2.8.5)

project(Geometry CXX)

add_library(Geometry
            Geometry.C
            DiscreteGeometry.C
)

removeWarnings( "Geometry.C;DiscreteGeometry.C")

INSTALL(TARGETS Geometry
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
