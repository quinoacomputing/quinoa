cmake_minimum_required(VERSION 2.8.5)

project(${QUINOA_EXECUTABLE} CXX)

# Basic configuration
add_definitions()

add_executable(${QUINOA_EXECUTABLE}
               Init.C
               TPLInfo/MKL.C
               TPLInfo/Boost.C
               TPLInfo/OpenMP.C
               TPLInfo/Zoltan.C
               TPLInfo/Silo.C
               TPLInfo/HDF5.C
               TPLInfo/Zlib.C
               RNGDriver.C
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
                      Paradigm
                      LinearAlgebra
                      Geometry
                      MonteCarlo
                      Base
                      IO
                      Mesh
                      Statistics
                      QuinoaControl
                      RNG
                      ${SILO_LIBRARY}
                      ${HDF5_LIBRARY}
                      ${Z_LIBRARY}
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_CORE_LIBRARY}
                      ${MKL_THREADED_LAYER_LIBRARY}
                      ${INTEL_OMP_RUNTIME_LIBRARY}
                      ${RNGSSE_LIBRARY}
                      ${BOOST_SYSTEM_LIBRARY}
)

INSTALL(TARGETS ${QUINOA_EXECUTABLE}
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

add_executable(${RNGTEST_EXECUTABLE}
               Init.C
               TPLInfo/MKL.C
               TPLInfo/Boost.C
               TPLInfo/OpenMP.C
               RNGDriver.C
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
                      #${TESTU01_MYLIB_LIBRARY}
                      ${TESTU01_PROBDIST_LIBRARY}
                      ${TESTU01_LIBRARY}
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_CORE_LIBRARY}
                      ${MKL_THREADED_LAYER_LIBRARY}
                      ${INTEL_OMP_RUNTIME_LIBRARY}
                      ${RNGSSE_LIBRARY}
                      ${BOOST_SYSTEM_LIBRARY}
                      #${TBB_LIBRARY}
)

INSTALL(TARGETS ${RNGTEST_EXECUTABLE}
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

add_executable(${GMSH2EXO_EXECUTABLE}
               Init.C
               TPLInfo/MKL.C
               TPLInfo/Boost.C
               TPLInfo/OpenMP.C
               Gmsh2ExoDriver.C
               Gmsh2Exo.C
)

if(MPI_COMPILE_FLAGS)
  set_target_properties(${GMSH2EXO_EXECUTABLE} PROPERTIES
                        COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${GMSH2EXO_EXECUTABLE} PROPERTIES
                        LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

target_link_libraries(${GMSH2EXO_EXECUTABLE}
                      Base
                      Gmsh2ExoControl
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_CORE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${BOOST_SYSTEM_LIBRARY}
)

INSTALL(TARGETS ${GMSH2EXO_EXECUTABLE}
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)
