cmake_minimum_required(VERSION 2.8.5)

project(${QUINOA_EXECUTABLE} CXX)

# Basic configuration
add_definitions()

add_executable(${QUINOA_EXECUTABLE}
               Init.C
               InitQuinoa.C
               Driver.C
               QuinoaDriver.C
               Quinoa.C
)

if(MPI_COMPILE_FLAGS)
  set_target_properties(${QUINOA_EXECUTABLE} PROPERTIES
                        COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${QUINOA_EXECUTABLE} PROPERTIES
                        LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

target_link_libraries(${QUINOA_EXECUTABLE}
                      Base
                      Paradigm
                      LinearAlgebra
                      Geometry
                      Physics
                      Model
                      IO
                      Mesh
                      Statistics
                      QuinoaControl
                      RNG
                      ${SILO_LIBRARY}
                      ${HDF5_LIBRARY}
                      ${Z_LIBRARY}
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_THREAD_LIBRARY}
                      ${MKL_CORE_LIBRARY}
                      ${INTEL_OMP_RUNTIME_LIBRARY}
                      #${MPI_LIBRARIES}
)

INSTALL(TARGETS ${QUINOA_EXECUTABLE}
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

add_executable(${RNGTEST_EXECUTABLE}
               Init.C
               InitRNGTest.C
               Driver.C
               RNGTestDriver.C
               RNGTest.C
)

if(MPI_COMPILE_FLAGS)
  set_target_properties(${RNGTEST_EXECUTABLE} PROPERTIES
                        COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${RNGTEST_EXECUTABLE} PROPERTIES
                        LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

target_link_libraries(${RNGTEST_EXECUTABLE}
                      Base
                      Paradigm
                      RNGTestControl
                      RNG
                      ${TESTU01_PROBDIST_LIBRARY}
                      #${TESTU01_MYLIB_LIBRARY}
                      ${TESTU01_LIBRARY}
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_THREAD_LIBRARY}
                      ${MKL_CORE_LIBRARY}
                      ${INTEL_OMP_RUNTIME_LIBRARY}
)

INSTALL(TARGETS ${RNGTEST_EXECUTABLE}
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)
