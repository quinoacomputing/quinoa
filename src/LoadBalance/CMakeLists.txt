project(LoadBalance CXX)

include(charm)

add_library(LoadBalance
            LinearMap.cpp
            UnsMeshMap.cpp
)

target_include_directories(LoadBalance PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${PROJECT_BINARY_DIR}/../LoadBalance
                           ${CHARM_INCLUDE_DIRS})

addCharmModule( "linearmap" "LoadBalance" )
addCharmModule( "unsmeshmap" "LoadBalance" )

set_target_properties(LoadBalance PROPERTIES
                      LIBRARY_OUTPUT_NAME quinoa_loadbalance)

INSTALL(TARGETS LoadBalance
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

add_library(ZoltanInterOp
            ZoltanInterOp.cpp)

target_include_directories(ZoltanInterOp PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${PROJECT_BINARY_DIR}/../Main
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_PATH}
                           ${Zoltan2_INCLUDE_DIRS}
                           ${BRIGAND_INCLUDE_DIRS})

set_target_properties(ZoltanInterOp PROPERTIES
                      LIBRARY_OUTPUT_NAME quinoa_zoltaninterop)

INSTALL(TARGETS ZoltanInterOp
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
