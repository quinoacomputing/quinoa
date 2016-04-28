ADD_LIBRARY(check_aec STATIC check_aec.c)
TARGET_LINK_LIBRARIES(check_aec aec)

ADD_EXECUTABLE(check_code_options check_code_options.c)
TARGET_LINK_LIBRARIES(check_code_options check_aec aec)

ADD_EXECUTABLE(check_buffer_sizes check_buffer_sizes.c)
TARGET_LINK_LIBRARIES(check_buffer_sizes check_aec aec)

ADD_EXECUTABLE(check_long_fs check_long_fs.c)
TARGET_LINK_LIBRARIES(check_long_fs check_aec aec)

ADD_TEST(check_code_options check_code_options)
ADD_TEST(check_buffer_sizes check_buffer_sizes)
ADD_TEST(check_long_fs check_long_fs)

ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND})
ADD_DEPENDENCIES(check
  check_code_options
  check_buffer_sizes
  check_long_fs)

IF(UNIX)
  ADD_CUSTOM_TARGET(bench
    COMMAND ${CMAKE_SOURCE_DIR}/tests/benc.sh
    COMMAND ${CMAKE_SOURCE_DIR}/tests/bdec.sh
    DEPENDS aec_client utime
    )
ENDIF(UNIX)
