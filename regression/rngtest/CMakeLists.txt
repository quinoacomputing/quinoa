# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

if (HAS_MKL)
  add_regression_test(SmallCrush_few_mkl ${RNGTEST_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES SmallCrush_few_mkl.q
                      ARGS -c SmallCrush_few_mkl.q -v)

  add_regression_test(Crush_mkl_mcg59 ${RNGTEST_EXECUTABLE}
                      NUMPES ${ManyPEs}
                      INPUTFILES Crush_mkl_mcg59.q
                      ARGS -c Crush_mkl_mcg59.q -v
                      LABELS stringent)

  add_regression_test(SmallCrush_all_mkl ${RNGTEST_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES SmallCrush_all_mkl.q
                      ARGS -c SmallCrush_all_mkl.q -v)
endif()


if (HAS_RNGSSE2)
  add_regression_test(SmallCrush_all_rngsse ${RNGTEST_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES SmallCrush_all_rngsse.q
                      ARGS -c SmallCrush_all_rngsse.q -v
                      LABELS stringent)

  add_regression_test(BigCrush_rngsse_mrg32k3a ${RNGTEST_EXECUTABLE}
                      NUMPES ${PROCESSOR_COUNT}
                      INPUTFILES BigCrush_rngsse_mrg32k3a.q
                      ARGS -c BigCrush_rngsse_mrg32k3a.q -v
                      LABELS stringent)
endif()

add_regression_test(SmallCrush_all_r123 ${RNGTEST_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES SmallCrush_all_r123.q
                    ARGS -c SmallCrush_all_r123.q -v)

add_regression_test(Crush_r123_threefry ${RNGTEST_EXECUTABLE}
                    NUMPES ${ManyPEs}
                    INPUTFILES Crush_r123_threefry.q
                    ARGS -c Crush_r123_threefry.q -v
                    LABELS stringent)

if (HAS_MKL AND HAS_RNGSSE2)
  add_regression_test(SmallCrush_mixed ${RNGTEST_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES SmallCrush_mixed.q
                      ARGS -c SmallCrush_mixed.q -v)
endif()
