### FileConv executable ########################################################

add_executable(${FILECONV_EXECUTABLE}
               FileConvDriver.cpp
               FileConv.cpp)

target_include_directories(${FILECONV_EXECUTABLE} PUBLIC
                           ${QUINOA_SOURCE_DIR}/IO
                           ${PROJECT_BINARY_DIR}/../Base)

config_executable(${FILECONV_EXECUTABLE})

target_link_libraries(${FILECONV_EXECUTABLE}
                      ${ROOTMESHIO}
                      ExodusIIMeshIO
                      Mesh
                      FileConvControl
                      Base
                      Config
                      Init
                      ${SEACASExodus_LIBRARIES}
                      ${ROOT_LIBRARIES}
                      ${NETCDF_LIBRARIES}       # only for static link
                      ${BACKWARD_LIBRARIES}
                      ${OMEGA_H_LIBRARIES}
                      ${LIBCXX_LIBRARIES}       # only for static link with libc++
                      ${LIBCXXABI_LIBRARIES})   # only for static link with libc++

# Add custom dependencies for FileConv's main Charm++ module
addCharmModule( "fileconv" "${FILECONV_EXECUTABLE}" )

add_dependencies( "fileconvCharmModule" "charestatecollectorCharmModule" )
