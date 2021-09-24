project(Base CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(Base
            Exception.cpp
            ProcessException.cpp
            ProcessControl.cpp
            LoadDistributor.cpp
            Timer.cpp
            Reader.cpp
            Writer.cpp
            PrintUtil.cpp
            ChareStateCollector.cpp
)

target_include_directories(Base PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${CHARM_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Main
                           ${BACKWARD_INCLUDE_DIRS}
                           ${BRIGAND_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_PATH}
                           ${PSTREAMS_INCLUDE_DIRS})

addCharmModule( "charestatecollector" "Base" )

set_target_properties(Base PROPERTIES LIBRARY_OUTPUT_NAME quinoa_base)

INSTALL(TARGETS Base
       RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
       LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
       ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
