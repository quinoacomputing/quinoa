# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(MassFracBeta ${WALKER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES massfracbeta.q stat.txt.std massfracbeta.ndiff.cfg
                    ARGS -c massfracbeta.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF massfracbeta.ndiff.cfg)

add_regression_test(MassFracBeta ${WALKER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES massfracbeta.q stat.txt.std massfracbeta.ndiff.cfg
                    ARGS -c massfracbeta.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF massfracbeta.ndiff.cfg)

add_regression_test(MassFracBeta ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES massfracbeta.q stat.txt.std massfracbeta.ndiff.cfg
                    ARGS -c massfracbeta.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF massfracbeta.ndiff.cfg)

add_regression_test(MassFracBeta_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES massfracbeta.q stat.txt.std massfracbeta.ndiff.cfg
                    ARGS -c massfracbeta.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF massfracbeta.ndiff.cfg)
