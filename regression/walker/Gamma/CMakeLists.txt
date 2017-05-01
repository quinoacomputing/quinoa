# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(Gamma ${WALKER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES gamma.q stat.txt.std gamma.ndiff.cfg
                    ARGS -c gamma.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF gamma.ndiff.cfg)

add_regression_test(Gamma ${WALKER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES gamma.q stat.txt.std gamma.ndiff.cfg
                    ARGS -c gamma.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF gamma.ndiff.cfg)

add_regression_test(Gamma ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gamma.q stat.txt.std gamma.ndiff.cfg
                    ARGS -c gamma.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF gamma.ndiff.cfg)

add_regression_test(Gamma_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gamma.q stat.txt.std gamma.ndiff.cfg
                    ARGS -c gamma.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF gamma.ndiff.cfg)
