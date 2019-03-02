### Walker executable ##########################################################

add_executable(${WALKER_EXECUTABLE}
               WalkerDriver.C
               WalkerPrint.C
               Walker.C)

config_executable(${WALKER_EXECUTABLE})

target_link_libraries(${WALKER_EXECUTABLE}
                      DiffEq
                      RNG
                      Walker
                      Statistics
                      IO
                      WalkerControl
                      Base
                      Config
                      Init
                      ${SEACASExodus_LIBRARIES}
                      ${LAPACKE_LIBRARIES}      # only if MKL not found
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${MKL_CORE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${RNGSSE2_LIBRARIES}
                      ${NETCDF_LIBRARIES}       # only for static link
                      ${HDF5_HL_LIBRARIES}      # only for static link
                      ${HDF5_C_LIBRARIES}
                      ${AEC_LIBRARIES}          # only for static link
                      ${BACKWARD_LIBRARIES})

# Add custom dependencies for Walker's main Charm++ module
addCharmModule( "walker" "${WALKER_EXECUTABLE}" )

add_dependencies( "walkerCharmModule" "charestatecollectorCharmModule" )
