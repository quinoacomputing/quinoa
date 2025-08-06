# Inciter executable ##########################################################

# set(INCITER_DRIVER_OBJ ${PROJECT_BINARY_DIR}/CMakeFiles/inciter.dir/InciterDriver.cpp.o)
# add_custom_command(
#     OUTPUT 
#         ${INCITER_DRIVER_OBJ} 
#     COMMAND
#         /projects/opt/centos8/x86_64/gcc/11.2.0/bin/g++
#         -DKOKKOS_DEPENDENCE 
#         -DQUINOA_CONFIG_MPI_ENABLED 
#         -I/vast/home/rjpark/quinoa/src 
#         -I/vast/home/rjpark/quinoa/src/Base 
#         -I/vast/home/rjpark/quinoa/src/Control 
#         -I/vast/home/rjpark/quinoa/src/Inciter 
#         -I/vast/home/rjpark/quinoa/src/Mesh 
#         -I/vast/home/rjpark/quinoa/src/PDE 
#         -I/vast/home/rjpark/quinoa-tpl/install/release/include/pegtl/include/tao 
#         -I/vast/home/rjpark/quinoa-tpl/install/release/charm/include 
#         -I/usr/include/lapacke -I/usr/include/cblas 
#         -I/vast/home/rjpark/quinoa/build/release/PDE/../Main 
#         -I/vast/home/rjpark/quinoa/src/IO 
#         -I/vast/home/rjpark/quinoa/src/Main 
#         -I/vast/home/rjpark/quinoa/src/LoadBalance 
#         -I/vast/home/rjpark/quinoa/src/Statistics 
#         -I/vast/home/rjpark/quinoa/src/LinearSolver 
#         -I/vast/home/rjpark/quinoa/src/Transfer 
#         -I/vast/home/rjpark/quinoa/build/release/Inciter/../Inciter 
#         -I/vast/home/rjpark/quinoa/build/release/Inciter/../LinearSolver 
#         -I/vast/home/rjpark/quinoa/build/release/Inciter/../Base 
#         -I/vast/home/rjpark/quinoa/build/release/Inciter/../IO 
#         -I/vast/home/rjpark/quinoa/build/release/Inciter/../Mesh 
#         -I/vast/home/rjpark/quinoa/build/release/Inciter/../Transfer 
#         -I/vast/home/rjpark/quinoa/build/release/Inciter/../Main 
#         -I/vast/home/rjpark/quinoa/build/release/Transfer 
#         -I/vast/home/rjpark/quinoa/build/release/Main 
#         -I/vast/home/rjpark/quinoa/build/release/IO 
#         -I/vast/home/rjpark/quinoa/src/Control/Inciter 
#         -I/vast/home/rjpark/quinoa/build/release/LoadBalance/../LoadBalance 
#         -I/vast/home/rjpark/quinoa/build/release/LoadBalance/../Main 
#         -I/vast/home/rjpark/quinoa/build/release/Base/../Main 
#         -I/vast/home/rjpark/quinoa/build/release/Main/../Main
#         -I/vast/home/rjpark/quinoa/build/release/IO/../Main 
#         -I/vast/home/rjpark/quinoa/build/release/IO/../IO
#         -I/vast/home/rjpark/quinoa/build/release/Mesh/../Main 
#         -I/vast/home/rjpark/quinoa/build/release/Statistics/../Main 
#         -I/vast/home/rjpark/quinoa/build/release/LinearSolver/../LinearSolver 
#         -isystem /vast/home/rjpark/quinoa-tpl/install/release/include 
#         -fdiagnostics-color -Wall -Wextra -Wcast-align -Wcast-qual 
#         -Wdisabled-optimization -Wfloat-equal 
#         -Wformat=2 -Wformat-nonliteral 
#         -Wformat-security -Wformat-y2k 
#         -Wimport -Winit-self -Winvalid-pch 
#         -Wmissing-field-initializers 
#         -Wmissing-format-attribute -Wmissing-noreturn 
#         -Wpacked -Wpointer-arith -Wredundant-decls 
#         -Wshadow -Wstack-protector -Wstrict-aliasing=2 
#         -Wunreachable-code -Wunused -Wunused-parameter 
#         -Wvariadic-macros -Wwrite-strings -Wno-sign-compare 
#         -Wno-unused-function -Wno-stack-protector 
#         -Wno-expansion-to-defined -Wno-int-in-bool-context
#         -Wno-cast-function-type -Wno-format-overflow 
#         -Wno-pragmas -Wno-unknown-pragmas 
#         -Wno-class-memaccess -O3 -DNDEBUG
#         -std=c++17 -extended-lambda 
#         #-Wext-lambda-captures-this 
#         #-arch=sm_61 
#         -MD -MT CMakeFiles/inciter.dir/InciterDriver.cpp.o 
#         -MF CMakeFiles/inciter.dir/InciterDriver.cpp.o.d 
#         -o CMakeFiles/inciter.dir/InciterDriver.cpp.o 
#         -c /vast/home/rjpark/quinoa/src/Main/InciterDriver.cpp
#     DEPENDS
#         charestatecollectorCharmModule
#         meshwriterCharmModule)

add_executable(${INCITER_EXECUTABLE}
               #${INCITER_DRIVER_OBJ}
               InciterDriver.cpp
               InciterPrint.cpp
               LBSwitch.cpp
               Inciter.cpp)

config_executable(${INCITER_EXECUTABLE})

message(STATUS "LUAAA ${LUA_LIBRARIES}")
target_link_libraries(${INCITER_EXECUTABLE}
                      InciterControl
                      PDE
                      Inciter
                      EOS
                      TransferDetails
                      TransportProblem
                      CGTransportPhysics
                      CompFlowProblem
                      MultiMatProblem
                      FVMultiMatPhysics
                      MultiSpeciesProblem
                      MultiSpeciesMixture
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
                      ${Trilinos_LIBRARIES}
                      ${LAPACKE_LIBRARIES}
                      ${CBLAS_LIBRARIES}
                      ${NETCDF_LIBRARIES}     # only for static link
                      ${BACKWARD_LIBRARIES}
                      ${LUA_LIBRARIES}
                      ${LIBCXX_LIBRARIES}     # only for static link with libc++
                      ${LIBCXXABI_LIBRARIES}) # only for static link with libc++

# Add custom dependencies for Inciter's main Charm++ module
addCharmModule( "inciter" "${INCITER_EXECUTABLE}" "-I${PROJECT_BINARY_DIR}")
addCharmModule( "lbswitch" "inciterCharmModule" )

add_dependencies( "inciterCharmModule" "charestatecollectorCharmModule" )
add_dependencies( "inciterCharmModule" "meshwriterCharmModule" )
