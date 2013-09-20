cmake_minimum_required(VERSION 2.8.5)

project(Physics CXX)

add_library(Physics
            Physics.C
            HomMix/HomMix.C
            HomHydro/HomHydro.C
            HomRT/HomRT.C
            SPINSFlow/SPINSFlow.C
)

INSTALL(TARGETS Physics
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
