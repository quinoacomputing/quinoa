cmake_minimum_required(VERSION 2.8.5)

project(Base CXX)

include(charm)

add_library(Base
            Exception.C
            ProcessException.C
            ProcessControl.C
            LoadDistributor.C
            Timer.C
            Reader.C
            Writer.C
            VectorReducer.C
            Table.C
)

target_include_directories(Base PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${CHARM_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Main
                           ${BACKWARDCPP_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_PATH}
                           ${PSTREAMS_INCLUDE_DIRS})

set_target_properties(Base PROPERTIES LIBRARY_OUTPUT_NAME quinoa_base)

INSTALL(TARGETS Base
       RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
       LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
       ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
