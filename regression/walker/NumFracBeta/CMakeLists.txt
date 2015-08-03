# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(NumFracBeta ${WALKER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES numfracbeta.q stat.txt.std numfracbeta.ndiff.cfg
                    ARGS -c numfracbeta.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF numfracbeta.ndiff.cfg)

add_regression_test(NumFracBeta ${WALKER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES numfracbeta.q stat.txt.std numfracbeta.ndiff.cfg
                    ARGS -c numfracbeta.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF numfracbeta.ndiff.cfg)

add_regression_test(NumFracBeta ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES numfracbeta.q stat.txt.std numfracbeta.ndiff.cfg
                    ARGS -c numfracbeta.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF numfracbeta.ndiff.cfg)

add_regression_test(NumFracBeta_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES numfracbeta.q stat.txt.std numfracbeta.ndiff.cfg
                    ARGS -c numfracbeta.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF numfracbeta.ndiff.cfg)

add_regression_test(NumFracBeta_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES numfracbeta.q stat.txt.std numfracbeta.ndiff.cfg
                    ARGS -c numfracbeta.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF numfracbeta.ndiff.cfg)
