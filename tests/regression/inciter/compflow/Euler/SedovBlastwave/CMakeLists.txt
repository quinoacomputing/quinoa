# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(compflow_euler_sedovblastwave_dgp1 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sedov_blastwave_dgp1.q unitsquare_01_3.6k.exo
                    ARGS -c sedov_blastwave_dgp1.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE sedov_blastwave_dgp1.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dgp1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sedov_blastwave_diag.ndiff.cfg)

# Parallel

add_regression_test(compflow_euler_sedovblastwave_dgp1_pe4 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES sedov_blastwave_dgp1.q unitsquare_01_3.6k.exo
                    ARGS -c sedov_blastwave_dgp1.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE sedov_blastwave_dgp1_pe4.std.exo.0
                                 sedov_blastwave_dgp1_pe4.std.exo.1
                                 sedov_blastwave_dgp1_pe4.std.exo.2
                                 sedov_blastwave_dgp1_pe4.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dgp1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sedov_blastwave_diag.ndiff.cfg)
