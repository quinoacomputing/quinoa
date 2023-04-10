## Inciter executable ##########################################################

add_executable(${INCITER_EXECUTABLE}
               InciterDriver.cpp
               InciterPrint.cpp
               LBSwitch.cpp
               Inciter.cpp)

config_executable(${INCITER_EXECUTABLE})

target_link_libraries(${INCITER_EXECUTABLE}
                      InciterControl
                      PDE
                      Inciter
                      EOS
                      TransportProblem
                      CGTransportPhysics
                      CompFlowProblem
                      MultiMatProblem
                      CGCompFlowPhysics
                      FVMultiMatPhysics
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
                      MeshDetect
                      Mesh
                      Statistics
                      LinearSolver
                      ${EXAM2M_LIBRARIES}
                      ${SEACASExodus_LIBRARIES}
                      ${Zoltan2_LIBRARIES}
                      ${LAPACKE_LIBRARIES}
                      ${NETCDF_LIBRARIES}     # only for static link
                      ${HDF5_HL_LIBRARIES}    # only for static link
                      ${HDF5_C_LIBRARIES}
                      ${BACKWARD_LIBRARIES}
                      ${LUA_LIBRARIES}
                      ${LIBCXX_LIBRARIES}     # only for static link with libc++
                      ${LIBCXXABI_LIBRARIES}) # only for static link with libc++

# Add custom dependencies for Inciter's main Charm++ module
addCharmModule( "inciter" "${INCITER_EXECUTABLE}" "-I${PROJECT_BINARY_DIR}")
addCharmModule( "lbswitch" "inciterCharmModule" )

add_dependencies( "inciterCharmModule" "charestatecollectorCharmModule" )
add_dependencies( "inciterCharmModule" "meshwriterCharmModule" )
