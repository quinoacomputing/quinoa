cmake_minimum_required(VERSION 2.8.5)

if(HAS_MKL)
  set(RNG "RNG" PARENT_SCOPE)

  project(RNG CXX)

  add_library(RNG
              MKLRNG.C)

  INSTALL(TARGETS RNG
          RUNTIME DESTINATION bin COMPONENT Runtime
          LIBRARY DESTINATION lib COMPONENT Runtime
          ARCHIVE DESTINATION lib COMPONENT Development)
endif()
