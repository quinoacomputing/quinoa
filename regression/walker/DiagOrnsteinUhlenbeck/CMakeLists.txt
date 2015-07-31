# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(DiagOrnsteinUhlenbeck ${WALKER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES diagou.q stat.txt.std diagou.ndiff.cfg
                    ARGS -c diagou.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

add_regression_test(DiagOrnsteinUhlenbeck ${WALKER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES diagou.q stat.txt.std diagou.ndiff.cfg
                    ARGS -c diagou.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

add_regression_test(DiagOrnsteinUhlenbeck ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES diagou.q stat.txt.std diagou.ndiff.cfg
                    ARGS -c diagou.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

add_regression_test(DiagOrnsteinUhlenbeck_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES diagou.q stat.txt.std diagou.ndiff.cfg
                    ARGS -c diagou.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

add_regression_test(DiagOrnsteinUhlenbeck_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES diagou.q stat.txt.std diagou.ndiff.cfg
                    ARGS -c diagou.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)
