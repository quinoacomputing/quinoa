cmake_minimum_required(VERSION 2.8.5)

project(${QUINOA_EXECUTABLE} CXX)

# Basic configuration
add_definitions()

# Include configuration/headers/libraries common to all code
INCLUDE(${QUINOA_SOURCE_DIR}/QuinoaCommon.cmake)

include_directories(${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Main
                    ${QUINOA_SOURCE_DIR}/Control
                    ${QUINOA_SOURCE_DIR}/Paradigm
                    ${QUINOA_SOURCE_DIR}/Model
                    ${QUINOA_SOURCE_DIR}/Physics
                    ${QUINOA_SOURCE_DIR}/Physics/HomogeneousDirichlet
                    ${QUINOA_SOURCE_DIR}/Physics/HomogeneousGeneralizedDirichlet
                    ${QUINOA_SOURCE_DIR}/Physics/SPINSFlow
)

add_executable(${QUINOA_EXECUTABLE}
               Driver.C
               Quinoa.C
)

target_link_libraries(${QUINOA_EXECUTABLE}
                      Base
                      LinearAlgebra
                      Mesh
                      Paradigm
                      HomogeneousDirichlet
                      IO
                      HomogeneousGeneralizedDirichlet
                      SPINSFlow
                      Statistics
                      Random
                      Dirichlet
                      GeneralizedDirichlet
                      MixModel
                      Model
                      Physics
                      ${MKL_LIBRARY}
                      ${MKL_THREAD_LIBRARY}
                      ${MKL_CORE_LIBRARY}
                      ${INTEL_OMP_LIBRARY}
                      ${PTHREAD_LIBRARY}
                      ${MATH_LIBRARY}
)

INSTALL(TARGETS ${QUINOA_EXECUTABLE}
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)
