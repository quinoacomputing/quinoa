### MeshConv executable ########################################################

add_executable(${MESHCONV_EXECUTABLE}
               MeshConvDriver.cpp
               MeshConv.cpp)

config_executable(${MESHCONV_EXECUTABLE})

target_include_directories(${MESHCONV_EXECUTABLE} PUBLIC
                           ${PROJECT_BINARY_DIR}/../Base)

target_link_libraries(${MESHCONV_EXECUTABLE}
                      NativeMeshIO
                      ExodusIIMeshIO
                      HyperMeshIO
                      MeshDetect
                      Mesh
                      MeshConvControl
                      Base
                      Config
                      Init
                      ${PUGIXML_LIBRARIES}
                      ${SEACASExodus_LIBRARIES}
                      ${NETCDF_LIBRARIES}       # only for static link
                      ${AEC_LIBRARIES}          # only for static link
                      ${BACKWARD_LIBRARIES}
                      ${OMEGA_H_LIBRARIES}
                      ${LIBCXX_LIBRARIES}       # only for static link with libc++
                      ${LIBCXXABI_LIBRARIES})   # only for static link with libc++

# Add custom dependencies for MeshConv's main Charm++ module
addCharmModule( "meshconv" "${MESHCONV_EXECUTABLE}" )

add_dependencies( "meshconvCharmModule" "charestatecollectorCharmModule" )
