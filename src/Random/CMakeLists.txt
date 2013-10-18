cmake_minimum_required(VERSION 2.8.5)

project(Random CXX)

add_library(Random
            MKL.C
            MKLRandom.C
            MKLRndStream.C
            MKLRndTable.C
)

INSTALL(TARGETS Random
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
