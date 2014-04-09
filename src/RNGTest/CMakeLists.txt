cmake_minimum_required(VERSION 2.8.5)

project(RNGTest CXX)

add_library(RNGTest
            TestU01Suite.C
            SmallCrush.C
            Crush.C
            BigCrush.C
)

INSTALL(TARGETS RNGTest
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
