cmake_minimum_required(VERSION 2.8.5)

project(RNG CXX)

if(HAS_MKL)
  set (MKLRNG_SOURCE "MKLRNG.C")
endif()

add_library(RNG
            TestU01Suite.C
            SmallCrush.C
            Crush.C
            BigCrush.C
            ${MKLRNG_SOURCE}
)

INSTALL(TARGETS RNG
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
