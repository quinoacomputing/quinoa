cmake_minimum_required(VERSION 2.8.5)

project(RNG CXX)

add_library(RNG
            RNGStack.C
)

target_include_directories(RNG PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${CHARM_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${BRIGAND_INCLUDE_DIRS}
                           ${TPL_INCLUDE_DIR}
                           ${RNGSSE2_INCLUDE_DIRS}
                           ${MKL_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Main)

set_target_properties(RNG PROPERTIES LIBRARY_OUTPUT_NAME quinoa_rng)

INSTALL(TARGETS RNG
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
