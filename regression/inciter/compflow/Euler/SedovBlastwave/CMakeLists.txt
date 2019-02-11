# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(compflow_euler_sedovblastwave_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sedov_blastwave_dg.q rectangle_01_1.5k.exo
                               diag_dg.std
                               exodiff_dg.cfg sedov_blastwave_dg.std.exo
                    ARGS -c sedov_blastwave_dg.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sedov_blastwave_dg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sedov_blastwave_diag.ndiff.cfg)
