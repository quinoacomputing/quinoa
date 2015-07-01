INCLUDE(TribitsAddExecutableAndTest)

include(${QT_USE_FILE})
  
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
if(UNIX AND (NOT APPLE) )
  message(STATUS "Target system is a non-apple, unix system. Testing ability to connect to X server")
  execute_process(
    COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/xtest.sh" 
    RESULT_VARIABLE CAN_CONNECT_TO_X
    OUTPUT_QUIET
    ERROR_QUIET)
  message(STATUS "Test Result: " ${CAN_CONNECT_TO_X})
  if(NOT ${CAN_CONNECT_TO_X} EQUAL 0)
    message(STATUS "Could not connect to an XServer. Disabling GUI unit tests by default.")
    SET(Optika_DO_GUI_UNIT_TESTS_DEFAULT "FALSE")
  else(NOT ${CAN_CONNECT_TO_X} EQUAL 0)
    message(STATUS "Found X server. Enabling GUI unit tests by default")
    SET(Optika_DO_GUI_UNIT_TESTS_DEFAULT "TRUE")
  endif(NOT ${CAN_CONNECT_TO_X} EQUAL 0)
else(UNIX AND (NOT APPLE) )
    message(STATUS "Not compling for a non-apple, unix system. Enabling GUI unit tests by default.")
    SET(Optika_DO_GUI_UNIT_TESTS_DEFAULT "TRUE")
endif(UNIX AND (NOT APPLE) )

SET(Optika_DO_GUI_UNIT_TESTS "${Optika_DO_GUI_UNIT_TESTS_DEFAULT}"
  CACHE BOOL "Build the Optika GUE unit tests.")

if(Optika_DO_GUI_UNIT_TESTS)

  message(STATUS "Enabling GUI unit tests.")

  QT4_GENERATE_MOC(GUI_UnitTests.cpp GUI_UnitTests.moc)
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
	  GUI_UnitTests
    SOURCES GUI_UnitTests.cpp GUI_UnitTests.moc
    FAIL_REGULAR_EXPRESSION "FAIL!"
    )
  
  TRIBITS_COPY_FILES_TO_BINARY_DIR(GUI_dep_test_copies
    DEST_FILES deptests.xml loadtest.xml loadtest.in.xml
    EXEDEPS GUI_UnitTests
  )

else(Optika_DO_GUI_UNIT_TESTS)

  message(STATUS "NOT enabling GUI unit tests.")

endif(Optika_DO_GUI_UNIT_TESTS)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  NonGUI_UnitTests
  SOURCES NonGUI_UnitTests.cpp
  ${TEUCHOS_STD_UNIT_TEST_MAIN}
)

