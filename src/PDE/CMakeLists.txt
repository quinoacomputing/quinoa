cmake_minimum_required(VERSION 2.8.5)

project(PDE CXX)

add_library(PDE
            PDEStack.C
)

INSTALL(TARGETS PDE
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
