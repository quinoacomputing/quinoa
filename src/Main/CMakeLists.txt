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
                    ${QUINOA_SOURCE_DIR}/Parser
                    ${QUINOA_SOURCE_DIR}/Random
                    ${QUINOA_SOURCE_DIR}/Model
                    ${QUINOA_SOURCE_DIR}/Model/Mass
                    ${QUINOA_SOURCE_DIR}/Model/Hydro
                    ${QUINOA_SOURCE_DIR}/Model/Mix
                    ${QUINOA_SOURCE_DIR}/Geometry
                    ${QUINOA_SOURCE_DIR}/Physics
                    ${QUINOA_SOURCE_DIR}/Physics/HomMix
                    ${QUINOA_SOURCE_DIR}/Physics/HomHydro
                    ${QUINOA_SOURCE_DIR}/Physics/HomRT
                    ${QUINOA_SOURCE_DIR}/Physics/SPINSFlow
)

add_executable(${QUINOA_EXECUTABLE}
               Driver.C
               Quinoa.C
)

target_link_libraries(${QUINOA_EXECUTABLE}
                      Base
                      Paradigm
                      LinearAlgebra
                      Geometry
                      HomMix
                      HomHydro
                      HomRT
                      SPINSFlow
                      IO
                      Mesh
                      Statistics
                      Physics
                      Beta
                      Dirichlet
                      GeneralizedDirichlet
                      SimplifiedLangevin
                      GeneralizedLangevin
                      Random
                      ${MKL_LIBRARY}
                      ${MKL_THREAD_LIBRARY}
                      ${MKL_CORE_LIBRARY}
                      ${INTEL_OMP_LIBRARY}
                      ${PTHREAD_LIBRARY}
                      ${MATH_LIBRARY}
                      Parser
)

INSTALL(TARGETS ${QUINOA_EXECUTABLE}
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)
