project(LinearSolver CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(LinearSolver
            CSR.cpp
            ConjugateGradients.cpp)

target_include_directories(LinearSolver PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/LinearSolver
                           ${CHARM_INCLUDE_DIRS}
                           ${BRIGAND_INCLUDE_DIRS}
                           ${HIGHWAYHASH_INCLUDE_DIRS}
                           ${BLAS_INCLUDE_DIRS}
                           ${LAPACKE_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../LinearSolver)

set_target_properties(LinearSolver PROPERTIES LIBRARY_OUTPUT_NAME
                      quinoa_linearsolver)

addCharmModule( "conjugategradients" "LinearSolver" )

INSTALL(TARGETS LinearSolver
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
