cmake_minimum_required(VERSION 2.8.5)

project(LinearSolver CXX)

add_library(LinearSolver
            CSR.cpp)

target_include_directories(LinearSolver PUBLIC
                           ${QUINOA_SOURCE_DIR}/Base)

set_target_properties(LinearSolver PROPERTIES LIBRARY_OUTPUT_NAME
                      quinoa_linearsolver)

INSTALL(TARGETS LinearSolver
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
