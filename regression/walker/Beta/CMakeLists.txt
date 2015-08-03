# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(Beta ${WALKER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES beta.q stat.txt.std beta.ndiff.cfg
                    ARGS -c beta.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF beta.ndiff.cfg)

add_regression_test(Beta ${WALKER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES beta.q stat.txt.std beta.ndiff.cfg
                    ARGS -c beta.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF beta.ndiff.cfg)

add_regression_test(Beta ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES beta.q stat.txt.std beta.ndiff.cfg
                    ARGS -c beta.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF beta.ndiff.cfg)

add_regression_test(Beta_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES beta.q stat.txt.std beta.ndiff.cfg
                    ARGS -c beta.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF beta.ndiff.cfg)

add_regression_test(Beta_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES beta.q stat.txt.std beta.ndiff.cfg
                    ARGS -c beta.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF beta.ndiff.cfg)
