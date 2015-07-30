# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(GeneralizedDirichlet ${WALKER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES gdir.q stat.txt.std gdir.ndiff.cfg
                    ARGS -c gdir.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF gdir.ndiff.cfg)

add_regression_test(GeneralizedDirichlet ${WALKER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES gdir.q stat.txt.std gdir.ndiff.cfg
                    ARGS -c gdir.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF gdir.ndiff.cfg)

add_regression_test(GeneralizedDirichlet ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gdir.q stat.txt.std gdir.ndiff.cfg
                    ARGS -c gdir.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF gdir.ndiff.cfg)

add_regression_test(GeneralizedDirichlet_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gdir.q stat.txt.std gdir.ndiff.cfg
                    ARGS -c gdir.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF gdir.ndiff.cfg)

add_regression_test(GeneralizedDirichlet_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES gdir.q stat.txt.std gdir.ndiff.cfg
                    ARGS -c gdir.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF gdir.ndiff.cfg)
