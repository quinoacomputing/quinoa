cmake_minimum_required(VERSION 2.8.5)

project(${QUINOA_EXECUTABLE} CXX)

# Basic configuration
add_definitions()

# Include configuration/headers/libraries common to all code
INCLUDE(${QUINOA_SOURCE_DIR}/QuinoaCommon.cmake)

include_directories(${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/IO
                    ${QUINOA_SOURCE_DIR}/LinearAlgebra
                    ${QUINOA_SOURCE_DIR}/Main
                    ${QUINOA_SOURCE_DIR}/Mesh
                    ${QUINOA_SOURCE_DIR}/Random
                    ${QUINOA_SOURCE_DIR}/Statistics
#                    ${QUINOA_SOURCE_DIR}/FEM
#                    ${QUINOA_SOURCE_DIR}/Particles
#                    ${QUINOA_SOURCE_DIR}/SDE
#                    ${QUINOA_SOURCE_DIR}/WallTreatment
)

add_executable(${QUINOA_EXECUTABLE}
               Quinoa.C
)

target_link_libraries(${QUINOA_EXECUTABLE}
                      Base
                      IO
                      LinearAlgebra
                      Mesh
                      Random
                      Statistics
#                      FEM
#                      Particles
#                      SDE
#                      WallTreatment
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
