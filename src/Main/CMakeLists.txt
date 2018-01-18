cmake_minimum_required(VERSION 2.8.5)

project(quinoa CXX)

### Quinoa executables build configuration #####################################

# Basic configuration
add_definitions()

# Include function for adding Charm++ modules
include(charm)

# Link executables with the charmc wrapper
STRING(REGEX REPLACE "<CMAKE_CXX_COMPILER>"
       "${CHARM_COMPILER} -module DistributedLB -c++ <CMAKE_CXX_COMPILER>"
       CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE}")

include(ConfigExecutable)

## UnitTest executable ########################################################
add_executable(${UNITTEST_EXECUTABLE}
               UnitTestDriver.C
               UnitTest.C
)

config_executable(${UNITTEST_EXECUTABLE})

target_link_libraries(${UNITTEST_EXECUTABLE}
                      Base
                      MeshIO
                      Mesh
                      RNG
                      LoadBalance
                      MeshRefinement
                      UnitTest
                      UnitTestControl
                      ${PUGIXML_LIBRARIES}
                      ${SEACASExodus_LIBRARIES}
                      ${RNGSSE2_LIBRARIES}
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${MKL_CORE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${Zoltan2_LIBRARIES}
                      ${H5PART_LIBRARIES}
                      ${ROOT_LIBRARIES}
                      ${NETCDF_LIBRARIES}       # only for static link
                      ${HDF5_HL_LIBRARIES}      # only for static link
                      ${HDF5_C_LIBRARIES}
                      ${AEC_LIBRARIES}          # only for static link
)

# Add custom dependencies for UnitTest's main Charm++ module
addCharmModule( "unittest" "${UNITTEST_EXECUTABLE}" )

### Inciter executable #########################################################
add_executable(${INCITER_EXECUTABLE}
               InciterDriver.C
               InciterPrint.C
               Inciter.C
)

config_executable(${INCITER_EXECUTABLE})

target_link_libraries(${INCITER_EXECUTABLE}
                      PDE
                      InciterControl
                      Inciter
                      LoadBalance
                      LinSys
                      Base
                      IO
                      MeshIO
                      Mesh
                      Statistics
                      MeshRefinement
                      ${PUGIXML_LIBRARIES}
                      ${SEACASExodus_LIBRARIES}
                      ${HYPRE_LIBRARIES}
                      ${ROOT_LIBRARIES}
                      ${Zoltan2_LIBRARIES}
                      ${RNGSSE2_LIBRARIES}
                      ${LAPACKE_LIBRARIES}      # only if MKL not found
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${MKL_CORE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${H5PART_LIBRARIES}
                      ${NETCDF_LIBRARIES}       # only for static link
                      ${HDF5_HL_LIBRARIES}      # only for static link
                      ${HDF5_C_LIBRARIES}
                      ${AEC_LIBRARIES}          # only for static link
)

# Add custom dependencies for Inciter's main Charm++ module
addCharmModule( "inciter" "${INCITER_EXECUTABLE}" )

### RNGTest executable #########################################################
if (HAS_TESTU01)

  add_executable(${RNGTEST_EXECUTABLE}
                 RNGTestDriver.C
                 RNGPrint.C
                 RNGTest.C
  )

  config_executable(${RNGTEST_EXECUTABLE})

  target_link_libraries(${RNGTEST_EXECUTABLE}
                        RNG
                        RNGTest
                        RNGTestControl
                        Base
                        ${TESTU01_LIBRARIES}
                        ${MKL_INTERFACE_LIBRARY}
                        ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                        ${MKL_CORE_LIBRARY}
                        ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                        ${RNGSSE2_LIBRARIES}
  )

  # Add custom dependencies for RNGTest's main Charm++ module
  addCharmModule( "rngtest" "${RNGTEST_EXECUTABLE}" )

endif(HAS_TESTU01)

### MeshConv executable ########################################################
add_executable(${MESHCONV_EXECUTABLE}
               MeshConvDriver.C
               MeshConv.C
)

config_executable(${MESHCONV_EXECUTABLE})

target_link_libraries(${MESHCONV_EXECUTABLE}
                      MeshIO
                      Mesh
                      MeshConvControl
                      Base
                      ${PUGIXML_LIBRARIES}
                      ${SEACASExodus_LIBRARIES}
                      ${H5PART_LIBRARIES}
                      ${ROOT_LIBRARIES}
                      ${NETCDF_LIBRARIES}       # only for static link
                      ${HDF5_HL_LIBRARIES}      # only for static link
                      ${HDF5_C_LIBRARIES}
                      ${AEC_LIBRARIES}          # only for static link
)

# Add custom dependencies for MeshConv's main Charm++ module
addCharmModule( "meshconv" "${MESHCONV_EXECUTABLE}" )

### Walker executable ##########################################################
add_executable(${WALKER_EXECUTABLE}
               WalkerDriver.C
               WalkerPrint.C
               Walker.C
)

config_executable(${WALKER_EXECUTABLE})

target_link_libraries(${WALKER_EXECUTABLE}
                      DiffEq
                      RNG
                      Walker
                      Statistics
                      IO
                      WalkerControl
                      Base
                      ${PUGIXML_LIBRARIES}
                      ${SEACASExodus_LIBRARIES}
                      ${LAPACKE_LIBRARIES}      # only if MKL not found
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${MKL_CORE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${RNGSSE2_LIBRARIES}
                      ${H5PART_LIBRARIES}
                      ${NETCDF_LIBRARIES}       # only for static link
                      ${HDF5_HL_LIBRARIES}      # only for static link
                      ${HDF5_C_LIBRARIES}
                      ${AEC_LIBRARIES}          # only for static link
)

# Add custom dependencies for Walker's main Charm++ module
addCharmModule( "walker" "${WALKER_EXECUTABLE}" )

### FileConv executable ########################################################
if (HAS_ROOT)

  add_executable(${FILECONV_EXECUTABLE}
                 FileConvDriver.C
                 FileConv.C
  )

  config_executable(${FILECONV_EXECUTABLE})

  target_link_libraries(${FILECONV_EXECUTABLE}
			MeshIO
			Mesh
                        FileConvControl
                        Base
			${PUGIXML_LIBRARIES}
	                ${SEACASExodus_LIBRARIES}
	                ${ROOT_LIBRARIES}
		        ${NETCDF_LIBRARIES}       # only for static link
  )

  # Add custom dependencies for FileConv's main Charm++ module
  addCharmModule( "fileconv" "${FILECONV_EXECUTABLE}" )

endif(HAS_ROOT)
