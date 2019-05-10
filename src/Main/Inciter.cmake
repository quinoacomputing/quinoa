## Inciter executable ##########################################################

add_executable(${INCITER_EXECUTABLE}
               InciterDriver.cpp
               InciterPrint.cpp
               LBSwitch.cpp
               Inciter.cpp)

config_executable(${INCITER_EXECUTABLE})

target_link_libraries(${INCITER_EXECUTABLE}
                      InciterControl
                      Inciter
                      PDE
                      TransportProblem
                      CGTransportPhysics
                      CompFlowProblem
                      CGCompFlowPhysics
                      Integrate
                      MeshRefinement
                      LoadBalance
                      ZoltanInterOp
                      Base
                      Config
                      Init
                      IO
                      MeshWriter
                      ExodusIIMeshIO
                      ${OMEGAHMESHIO}
                      ${ROOTMESHIO}
                      MeshDetect
                      Mesh
                      Statistics
                      ${SEACASExodus_LIBRARIES}
                      ${ROOT_LIBRARIES}
                      ${Zoltan2_LIBRARIES}
                      ${RNGSSE2_LIBRARIES}
                      ${LAPACKE_LIBRARIES}      # only if MKL not found
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${MKL_CORE_LIBRARY}
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${NETCDF_LIBRARIES}       # only for static link
                      ${HDF5_HL_LIBRARIES}      # only for static link
                      ${HDF5_C_LIBRARIES}
                      ${AEC_LIBRARIES}          # only for static link
                      ${BACKWARD_LIBRARIES}
                      ${OMEGA_H_LIBRARIES})

# Add custom dependencies for Inciter's main Charm++ module
addCharmModule( "inciter" "${INCITER_EXECUTABLE}" )
addCharmModule( "lbswitch" "${INCITER_EXECUTABLE}" )

add_dependencies( "inciterCharmModule" "charestatecollectorCharmModule" )
add_dependencies( "inciterCharmModule" "meshwriterCharmModule" )
