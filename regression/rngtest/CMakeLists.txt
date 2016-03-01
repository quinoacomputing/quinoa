# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

if (HAS_MKL)
  add_regression_test(SmallCrush_few_mkl ${RNGTEST_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES SmallCrush_few_mkl.q
                      ARGS -c SmallCrush_few_mkl.q -v)

  add_regression_test(Crush_MKLmcg59 ${RNGTEST_EXECUTABLE}
                      NUMPES ${ManyPEs}
                      INPUTFILES Crush_MKLmcg59.q
                      ARGS -c Crush_MKLmcg59.q -v
                      LABELS stringent)

  add_regression_test(SmallCrush_all_mkl ${RNGTEST_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES SmallCrush_all_mkl.q
                      ARGS -c SmallCrush_all_mkl.q -v)
endif()


add_regression_test(SmallCrush_all_rngsse ${RNGTEST_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES SmallCrush_all_rngsse.q
                    ARGS -c SmallCrush_all_rngsse.q -v
                    LABELS stringent)

add_regression_test(BigCrush_RNGSSEmrg32k3a ${RNGTEST_EXECUTABLE}
                    NUMPES ${PROCESSOR_COUNT}
                    INPUTFILES BigCrush_RNGSSEmrg32k3a.q
                    ARGS -c BigCrush_RNGSSEmrg32k3a.q -v
                    LABELS stringent)
