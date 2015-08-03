# Cmake code common to all regression testing

# Quinoa cmake modules directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# Save paths required for running tests
set(NDIFF_EXECUTABLE ${CMAKE_BINARY_DIR}/regression/ndiff/numdiff)
set(EXODIFF_EXECUTABLE ${TPL_DIR}/bin/exodiff)
set(TEST_RUNNER ${CMAKE_SOURCE_DIR}/../cmake/test_runner.cmake)

# Include function used to add regression tests
include(add_regression_test)

# Include regression test directories

# Walker regression tests
add_subdirectory(walker/Dirichlet)
add_subdirectory(walker/GeneralizedDirichlet)
add_subdirectory(walker/SkewNormal)
add_subdirectory(walker/DiagOrnsteinUhlenbeck)
add_subdirectory(walker/OrnsteinUhlenbeck)
add_subdirectory(walker/Beta)
add_subdirectory(walker/Gamma)
