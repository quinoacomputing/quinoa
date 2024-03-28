# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(compflow_euler_taylor_green_dgp2 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES taylor_green_dgp2.q unitcube_1k.exo
                    ARGS -c taylor_green_dgp2.q -i unitcube_1k.exo -v
                    BIN_BASELINE taylor_green_dgp2.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dgp2.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF taylor_green_diag.ndiff.cfg
                    LABELS dg)

add_regression_test(compflow_euler_taylor_green_dgp2_cfl ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES taylor_green_dgp2_cfl.q unitcube_1k.exo
                    ARGS -c taylor_green_dgp2_cfl.q -i unitcube_1k.exo -v
                    BIN_BASELINE taylor_green_dgp2_cfl.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dgp2_cfl.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF taylor_green_diag.ndiff.cfg
                    LABELS dg)

# Parallel + no virtualization

add_regression_test(compflow_euler_taylor_green_dgp2 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES taylor_green_dgp2.q unitcube_1k.exo
                    ARGS -c taylor_green_dgp2.q -i unitcube_1k.exo -v
                    BIN_BASELINE taylor_green_dgp2_pe4.std.exo.0
                                 taylor_green_dgp2_pe4.std.exo.1
                                 taylor_green_dgp2_pe4.std.exo.2
                                 taylor_green_dgp2_pe4.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    TEXT_BASELINE diag_dgp2.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF taylor_green_diag.ndiff.cfg
                    LABELS dg)

add_regression_test(compflow_euler_taylor_green_dgp2_cfl ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES taylor_green_dgp2_cfl.q unitcube_1k.exo
                    ARGS -c taylor_green_dgp2_cfl.q -i unitcube_1k.exo -v
                    BIN_BASELINE taylor_green_dgp2_cfl_pe4.std.exo.0
                                 taylor_green_dgp2_cfl_pe4.std.exo.1
                                 taylor_green_dgp2_cfl_pe4.std.exo.2
                                 taylor_green_dgp2_cfl_pe4.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    TEXT_BASELINE diag_dgp2_cfl.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF taylor_green_diag.ndiff.cfg
                    LABELS dg)

# Parallel + virtualization
