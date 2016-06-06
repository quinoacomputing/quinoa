cmake_minimum_required(VERSION 2.4)

PROJECT(numdiff C)

if(NOT "${CMAKE_C_FLAGS}" MATCHES "std=c99")
   if(NOT (WIN32 AND CMAKE_C_COMPILER_ID STREQUAL "Intel"))
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
   else()
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Qstd=c99")
   endif()
endif()

file(GLOB numdiff_src src/*.c)

include_directories(src)
add_executable(${NUMDIFF_EXECUTABLE} ${numdiff_src})

# Conditionally enforce static linking
if(NOT BUILD_SHARED_LIBS)
  set_target_properties(${NUMDIFF_EXECUTABLE} PROPERTIES LINK_FLAGS "-static")
endif()

if(UNIX)
   target_link_libraries(${NUMDIFF_EXECUTABLE} m)
endif()

install(TARGETS ${NUMDIFF_EXECUTABLE}
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development)
