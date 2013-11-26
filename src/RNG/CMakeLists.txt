cmake_minimum_required(VERSION 2.8.5)

project(RNG CXX)

if (HAS_MKL)
  set (MKL_SOURCE "MKL.C")
  set (MKLRNG_SOURCE "MKLRNG.C")
  set (MKLRNGWRAPPERS_SOURCE "MKLRNGWrappers.C")
  set (MKLRNDSTREAM_SOURCE "MKLRndStream.C")
  set (MKLRNDTABLE_SOURCE "MKLRndTable.C")
endif()

add_library(RNG
            RNG.C
            SmallCrush.C
            TestU01Suite.C
            ${MKL_SOURCE}
            ${MKLRNG_SOURCE}
            ${MKLRNGWRAPPERS_SOURCE}
            ${MKLRNDSTREAM_SOURCE}
            ${MKLRNDTABLE_SOURCE}
)

INSTALL(TARGETS RNG
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
