# See regression/CMakeLists.txt for documentation on the arguments to
# add_regression_test().

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
