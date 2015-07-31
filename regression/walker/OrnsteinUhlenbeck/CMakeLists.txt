# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(OrnsteinUhlenbeck ${WALKER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES ou.q stat.txt.std ou.ndiff.cfg
                    ARGS -c ou.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF ou.ndiff.cfg)

add_regression_test(OrnsteinUhlenbeck ${WALKER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES ou.q stat.txt.std ou.ndiff.cfg
                    ARGS -c ou.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF ou.ndiff.cfg)

add_regression_test(OrnsteinUhlenbeck ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES ou.q stat.txt.std ou.ndiff.cfg
                    ARGS -c ou.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF ou.ndiff.cfg)

add_regression_test(OrnsteinUhlenbeck_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES ou.q stat.txt.std ou.ndiff.cfg
                    ARGS -c ou.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF ou.ndiff.cfg)

add_regression_test(OrnsteinUhlenbeck_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES ou.q stat.txt.std ou.ndiff.cfg
                    ARGS -c ou.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF ou.ndiff.cfg)
