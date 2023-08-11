# Include function for adding Charm++ modules
include(charm)

add_library(Worker
            Worker.cpp
            M2MTransfer.cpp)

target_include_directories(Worker PUBLIC
                           ${PROJECT_SOURCE_DIR}
                           ${PROJECT_SOURCE_DIR}/Base
                           ${PROJECT_SOURCE_DIR}/Mesh
                           ${PROJECT_BINARY_DIR}/Transfer
                           ${PROJECT_BINARY_DIR}/Main
                           ${PROJECT_BINARY_DIR}/IO
                           ${QUINOA_SOURCE_DIR}/Control
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${HIGHWAYHASH_INCLUDE_DIRS}
)

set_target_properties(Worker PROPERTIES LIBRARY_OUTPUT_NAME exam2m_worker)

INSTALL(TARGETS Worker
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

#INSTALL(DIRECTORY "${PROJECT_SOURCE_DIR}/Transfer/"
#        DESTINATION "include"
#        FILES_MATCHING PATTERN "M2MTransfer.hpp")
#INSTALL(DIRECTORY "${PROJECT_SOURCE_DIR}/NoWarning/"
#        DESTINATION "include/NoWarning/"
#        FILES_MATCHING PATTERN "m2mtransfer.de*.h")
#INSTALL(DIRECTORY "${PROJECT_SOURCE_DIR}/NoWarning/"
#        DESTINATION "include/NoWarning/"
#        FILES_MATCHING PATTERN "worker.de*.h")
#
#INSTALL(DIRECTORY "${PROJECT_BINARY_DIR}/Transfer/"
#        DESTINATION "include/Transfer/"
#        FILES_MATCHING PATTERN "m2mtransfer.de*.h")
#INSTALL(DIRECTORY "${PROJECT_BINARY_DIR}/Transfer/"
#        DESTINATION "include/Transfer/"
#        FILES_MATCHING PATTERN "worker.de*.h")

addCharmModule( "worker" "Worker" )
addCharmModule( "m2mtransfer" "Worker" )
