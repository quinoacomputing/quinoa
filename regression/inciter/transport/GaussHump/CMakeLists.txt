# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Virtualization = 0

add_regression_test(gauss_hump ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES gauss_hump.q unitsquare_01_3.6k.exo diag.std
                               exodiff.cfg gauss_hump.std.exo
                    ARGS -c gauss_hump.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE gauss_hump.std.exo
                    BIN_RESULT out.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)

add_regression_test(gauss_hump ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gauss_hump.q unitsquare_01_3.6k.exo diag.std
                               exodiff.cfg
                               gauss_hump_pe4_u0.0.std.exo.0
                               gauss_hump_pe4_u0.0.std.exo.1
                               gauss_hump_pe4_u0.0.std.exo.2
                               gauss_hump_pe4_u0.0.std.exo.3
                    ARGS -c gauss_hump.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE gauss_hump_pe4_u0.0.std.exo.0
                                 gauss_hump_pe4_u0.0.std.exo.1
                                 gauss_hump_pe4_u0.0.std.exo.2
                                 gauss_hump_pe4_u0.0.std.exo.3
                    BIN_RESULT out.0
                               out.1
                               out.2
                               out.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)

# Parallel + virtualization

add_regression_test(gauss_hump_u0.4 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gauss_hump.q unitsquare_01_3.6k.exo diag.std
                               exodiff.cfg
                               gauss_hump_pe4_u0.5.std.exo.0
                               gauss_hump_pe4_u0.5.std.exo.1
                               gauss_hump_pe4_u0.5.std.exo.2
                               gauss_hump_pe4_u0.5.std.exo.3
                               gauss_hump_pe4_u0.5.std.exo.4
                               gauss_hump_pe4_u0.5.std.exo.5
                               gauss_hump_pe4_u0.5.std.exo.6
                               gauss_hump_pe4_u0.5.std.exo.7
                    ARGS -c gauss_hump.q -i unitsquare_01_3.6k.exo -v -u 0.5
                    BIN_BASELINE gauss_hump_pe4_u0.5.std.exo.0
                                 gauss_hump_pe4_u0.5.std.exo.1
                                 gauss_hump_pe4_u0.5.std.exo.2
                                 gauss_hump_pe4_u0.5.std.exo.3
                                 gauss_hump_pe4_u0.5.std.exo.4
                                 gauss_hump_pe4_u0.5.std.exo.5
                                 gauss_hump_pe4_u0.5.std.exo.6
                                 gauss_hump_pe4_u0.5.std.exo.7
                    BIN_RESULT out.0
                               out.1
                               out.2
                               out.3
                               out.4
                               out.5
                               out.6
                               out.7
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)
