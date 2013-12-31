cmake_minimum_required(VERSION 2.8.5)

project(MonteCarlo CXX)

add_library(MonteCarlo
            TestSDE.C
            Physics.C
            HomMix.C
            HomHydro.C
            HomRT.C
            SPINSFlow.C
)

INSTALL(TARGETS MonteCarlo
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
