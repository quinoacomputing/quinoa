cmake_minimum_required(VERSION 2.8.5)

project(Statistics CXX)

add_library(Statistics
            Statistics.C
)

removeWarnings( "Statistics.C" )

INSTALL(TARGETS Statistics
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
