cmake_minimum_required(VERSION 2.8.5)

project(SDE CXX)

add_library(SDE
            Dirichlet.C
            GenDirichlet.C
)

INSTALL(TARGETS SDE
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
