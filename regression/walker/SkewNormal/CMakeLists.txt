# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(SkewNormal ${WALKER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES skew.q stat.txt.std skew.ndiff.cfg
                    ARGS -c skew.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

add_regression_test(SkewNormal ${WALKER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES skew.q stat.txt.std skew.ndiff.cfg
                    ARGS -c skew.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

add_regression_test(SkewNormal ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES skew.q stat.txt.std skew.ndiff.cfg
                    ARGS -c skew.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

add_regression_test(SkewNormal_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES skew.q stat.txt.std skew.ndiff.cfg
                    ARGS -c skew.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

add_regression_test(SkewNormal_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES skew.q stat.txt.std skew.ndiff.cfg
                    ARGS -c skew.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF skew.ndiff.cfg)
