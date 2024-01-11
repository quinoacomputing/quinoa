# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(compflow_euler_timedepbc ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES timedep_bc.q rectangle_01_1.5k.exo
                    ARGS -c timedep_bc.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE timedep_bc.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF timedep_bc.ndiff.cfg
                    LABELS alecg)

add_regression_test(compflow_euler_timedepbc_oversetfe ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES timedep_bc_oversetfe.q rectangle_01_1.5k.exo
                    ARGS -c timedep_bc_oversetfe.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE timedep_bc.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF timedep_bc.ndiff.cfg
                    LABELS oversetfe)
