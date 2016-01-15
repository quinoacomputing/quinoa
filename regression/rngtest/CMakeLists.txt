# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(SmallCrush_few ${RNGTEST_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES SmallCrush_few.q
                    ARGS -c SmallCrush_few.q -v)

add_regression_test(SmallCrush_all ${RNGTEST_EXECUTABLE}
                    NUMPES ${ManyPEs}
                    INPUTFILES SmallCrush_all.q
                    ARGS -c SmallCrush_all.q -v)

add_regression_test(Crush_MKLmcg59 ${RNGTEST_EXECUTABLE}
                    NUMPES ${ManyPEs}
                    INPUTFILES Crush_MKLmcg59.q
                    ARGS -c Crush_MKLmcg59.q -v)

add_regression_test(BigCrush_RNGSSEmrg32k3a ${RNGTEST_EXECUTABLE}
                    NUMPES ${PROCESSOR_COUNT}
                    INPUTFILES BigCrush_RNGSSEmrg32k3a.q
                    ARGS -c BigCrush_RNGSSEmrg32k3a.q -v)
