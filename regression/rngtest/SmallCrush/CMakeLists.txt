# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(few ${RNGTEST_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES few.q
                    ARGS -c few.q -v)

add_regression_test(all ${RNGTEST_EXECUTABLE}
                    NUMPES ${PROCESSOR_COUNT}
                    INPUTFILES all.q
                    ARGS -c all.q -v)
