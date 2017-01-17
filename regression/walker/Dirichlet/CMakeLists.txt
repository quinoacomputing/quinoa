# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

if (HAS_MKL)
  add_regression_test(Dirichlet_mkl ${WALKER_EXECUTABLE}
                      NUMPES 1
                      INPUTFILES dir_mkl.q stat.txt.std dir.ndiff.cfg
                      ARGS -c dir_mkl.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF dir.ndiff.cfg)

  add_regression_test(Dirichlet_mkl ${WALKER_EXECUTABLE}
                      NUMPES 2
                      INPUTFILES dir_mkl.q stat.txt.std dir.ndiff.cfg
                      ARGS -c dir_mkl.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF dir.ndiff.cfg)

  add_regression_test(Dirichlet_mkl ${WALKER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES dir_mkl.q stat.txt.std dir.ndiff.cfg
                      ARGS -c dir_mkl.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF dir.ndiff.cfg)

  add_regression_test(Dirichlet_mkl_u0.9 ${WALKER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES dir_mkl.q stat.txt.std dir.ndiff.cfg
                      ARGS -c dir_mkl.q -v -u 0.9
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF dir.ndiff.cfg)

  add_regression_test(Dirichlet_mkl_u0.9 ${WALKER_EXECUTABLE}
                      NUMPES 8
                      INPUTFILES dir_mkl.q stat.txt.std dir.ndiff.cfg
                      ARGS -c dir_mkl.q -v -u 0.9
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF dir.ndiff.cfg)
endif()

if (HAS_RNGSSE2)
  add_regression_test(Dirichlet_rngsse ${WALKER_EXECUTABLE}
                      NUMPES 1
                      INPUTFILES dir_rngsse.q stat.txt.std dir.ndiff.cfg
                      ARGS -c dir_rngsse.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF dir.ndiff.cfg)

  add_regression_test(Dirichlet_rngsse ${WALKER_EXECUTABLE}
                      NUMPES 2
                      INPUTFILES dir_rngsse.q stat.txt.std dir.ndiff.cfg
                      ARGS -c dir_rngsse.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF dir.ndiff.cfg)

  add_regression_test(Dirichlet_rngsse ${WALKER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES dir_rngsse.q stat.txt.std dir.ndiff.cfg
                      ARGS -c dir_rngsse.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF dir.ndiff.cfg)

  add_regression_test(Dirichlet_rngsse_u0.9 ${WALKER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES dir_rngsse.q stat.txt.std dir.ndiff.cfg
                      ARGS -c dir_rngsse.q -v -u 0.9
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF dir.ndiff.cfg)

  add_regression_test(Dirichlet_rngsse_u0.9 ${WALKER_EXECUTABLE}
                      NUMPES 8
                      INPUTFILES dir_rngsse.q stat.txt.std dir.ndiff.cfg
                      ARGS -c dir_rngsse.q -v -u 0.9
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF dir.ndiff.cfg)
endif()

add_regression_test(Dirichlet ${WALKER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES dir.q stat.txt.std dir.ndiff.cfg
                    ARGS -c dir.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF dir.ndiff.cfg)

add_regression_test(Dirichlet ${WALKER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES dir.q stat.txt.std dir.ndiff.cfg
                    ARGS -c dir.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF dir.ndiff.cfg)

add_regression_test(Dirichlet ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES dir.q stat.txt.std dir.ndiff.cfg
                    ARGS -c dir.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF dir.ndiff.cfg)

add_regression_test(Dirichlet_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES dir.q stat.txt.std dir.ndiff.cfg
                    ARGS -c dir.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF dir.ndiff.cfg)

add_regression_test(Dirichlet_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES dir.q stat.txt.std dir.ndiff.cfg
                    ARGS -c dir.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF dir.ndiff.cfg)
