cmake_minimum_required(VERSION 2.8.5)

project(RNG CXX)

add_library(RNG
            RNGStack.C
)

INSTALL(TARGETS RNG
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
