# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(MixMassFracBeta ${WALKER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES mixmassfracbeta.q stat.txt.std
                               mixmassfracbeta.ndiff.cfg
                    ARGS -c mixmassfracbeta.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF mixmassfracbeta.ndiff.cfg)

add_regression_test(MixMassFracBeta ${WALKER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES mixmassfracbeta.q stat.txt.std
                               mixmassfracbeta.ndiff.cfg
                    ARGS -c mixmassfracbeta.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF mixmassfracbeta.ndiff.cfg)

add_regression_test(MixMassFracBeta ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES mixmassfracbeta.q stat.txt.std
                               mixmassfracbeta.ndiff.cfg
                    ARGS -c mixmassfracbeta.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF mixmassfracbeta.ndiff.cfg)

add_regression_test(MixMassFracBeta_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES mixmassfracbeta.q stat.txt.std
                               mixmassfracbeta.ndiff.cfg
                    ARGS -c mixmassfracbeta.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF mixmassfracbeta.ndiff.cfg)

add_regression_test(MixMassFracBeta_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES mixmassfracbeta.q stat.txt.std
                               mixmassfracbeta.ndiff.cfg
                    ARGS -c mixmassfracbeta.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF mixmassfracbeta.ndiff.cfg)
