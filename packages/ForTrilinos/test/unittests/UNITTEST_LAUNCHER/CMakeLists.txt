
INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  )

TRIBITS_ADD_EXECUTABLE(unittest_launcher_serial SOURCES unittest_launcher_serial.cpp)

SET_TARGET_PROPERTIES(
  ForTrilinos_unittest_launcher_serial PROPERTIES
  LINKER_LANGUAGE CXX
  )

