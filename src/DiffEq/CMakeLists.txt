cmake_minimum_required(VERSION 2.8.5)

project(DiffEq CXX)

add_library(DiffEq
            DiffEqStack.C
)

INSTALL(TARGETS DiffEq
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
