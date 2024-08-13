# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(compflow_euler_nleg_alecg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES nleg_alecg.lua unitcube_1k.exo
                    ARGS -c nleg_alecg.lua -i unitcube_1k.exo -v
                    BIN_BASELINE nleg_alecg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_alecg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF nleg_diag.ndiff.cfg
                    LABELS alecg)

add_regression_test(compflow_euler_nleg_dgp2 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES nleg_dgp2.lua unitcube_1k.exo
                    ARGS -c nleg_dgp2.lua -i unitcube_1k.exo -v
                    BIN_BASELINE nleg_dgp2.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dgp2.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF nleg_diag.ndiff.cfg
                    LABELS dg)

# Parallel + no virtualization

add_regression_test(compflow_euler_nleg_dgp2 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES nleg_dgp2.lua unitcube_1k.exo
                    ARGS -c nleg_dgp2.lua -i unitcube_1k.exo -v
                    BIN_BASELINE nleg_dgp2_pe4.std.exo.0
                                 nleg_dgp2_pe4.std.exo.1
                                 nleg_dgp2_pe4.std.exo.2
                                 nleg_dgp2_pe4.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    TEXT_BASELINE diag_dgp2.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF nleg_diag.ndiff.cfg
                    LABELS dg)

# Parallel + virtualization

add_regression_test(compflow_euler_nleg_alecg_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES nleg_alecg.lua unitcube_1k.exo
                    ARGS -c nleg_alecg.lua -i unitcube_1k.exo -v -u 0.5
                    BIN_BASELINE nleg_alecg_pe4_u0.5.std.exo.0
                                 nleg_alecg_pe4_u0.5.std.exo.1
                                 nleg_alecg_pe4_u0.5.std.exo.2
                                 nleg_alecg_pe4_u0.5.std.exo.3
                                 nleg_alecg_pe4_u0.5.std.exo.4
                                 nleg_alecg_pe4_u0.5.std.exo.5
                                 nleg_alecg_pe4_u0.5.std.exo.6
                                 nleg_alecg_pe4_u0.5.std.exo.7
                    BIN_RESULT out.e-s.0.8.0
                               out.e-s.0.8.1
                               out.e-s.0.8.2
                               out.e-s.0.8.3
                               out.e-s.0.8.4
                               out.e-s.0.8.5
                               out.e-s.0.8.6
                               out.e-s.0.8.7
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_alecg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF nleg_diag.ndiff.cfg
                    LABELS alecg)
