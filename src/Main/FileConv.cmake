### FileConv executable ########################################################

add_executable(${FILECONV_EXECUTABLE}
               FileConvDriver.C
               FileConv.C)

target_include_directories(${FILECONV_EXECUTABLE} PUBLIC
                           ${QUINOA_SOURCE_DIR}/IO
                           ${PROJECT_BINARY_DIR}/../Base)

config_executable(${FILECONV_EXECUTABLE})

target_link_libraries(${FILECONV_EXECUTABLE}
                      ExodusIIMeshIO
                      ${ROOTMESHIO}
                      Mesh
                      FileConvControl
                      Base
                      Config
                      Init
                      ${SEACASExodus_LIBRARIES}
                      ${ROOT_LIBRARIES}
                      ${NETCDF_LIBRARIES}       # only for static link
                      ${BACKWARD_LIBRARIES}
                      ${OMEGA_H_LIBRARIES})

# Add custom dependencies for FileConv's main Charm++ module
addCharmModule( "fileconv" "${FILECONV_EXECUTABLE}" )

add_dependencies( "fileconvCharmModule" "charestatecollectorCharmModule" )
