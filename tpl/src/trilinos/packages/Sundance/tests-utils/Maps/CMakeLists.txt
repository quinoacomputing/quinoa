# CMake tests specification 


INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
        mapPerfTest
        SOURCES mapPerfTest.cpp
        COMM serial
)
