### RNGTest executable #########################################################

add_executable(${RNGTEST_EXECUTABLE}
               RNGTestDriver.cpp
               RNGPrint.cpp
               RNGTest.cpp)

target_include_directories(${RNGTEST_EXECUTABLE} PUBLIC
                           ${QUINOA_SOURCE_DIR}/RNGTest
                           ${QUINOA_SOURCE_DIR}/Main)

config_executable(${RNGTEST_EXECUTABLE})

target_link_libraries(${RNGTEST_EXECUTABLE}
                      RNG
                      RNGTest
                      RNGTestControl
                      Base
                      Config
                      Init
                      ${TESTU01_LIBRARIES}
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${MKL_CORE_LIBRARY}
                      ${MKL_INTERFACE_LIBRARY}
                      ${MKL_SEQUENTIAL_LAYER_LIBRARY}
                      ${RNGSSE2_LIBRARIES}
                      ${BACKWARD_LIBRARIES})

# Add custom dependencies for RNGTest's main Charm++ module
addCharmModule( "rngtest" "${RNGTEST_EXECUTABLE}" )

add_dependencies( "rngtestCharmModule" "charestatecollectorCharmModule" )
