# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Virtualization = 0

add_regression_test(gauss_hump_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES gauss_hump.q unitsquare_01_3.6k.exo diag.std
                               exodiff.cfg gauss_hump.std.exo
                    ARGS -c gauss_hump.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE gauss_hump.std.exo
                    BIN_RESULT out.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)

add_regression_test(gauss_hump_dg ${INCITER_EXECUTABLE}
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

add_regression_test(gauss_hump_dg_u0.5 ${INCITER_EXECUTABLE}
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

# Parallel + virtualization + migration

add_regression_test(gauss_hump_dg_u0.9_migr ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gauss_hump.q unitsquare_01_3.6k.exo diag.std
                               exodiff.cfg
                               gauss_hump_pe4_u0.9.std.exo.0
                               gauss_hump_pe4_u0.9.std.exo.1
                               gauss_hump_pe4_u0.9.std.exo.2
                               gauss_hump_pe4_u0.9.std.exo.3
                               gauss_hump_pe4_u0.9.std.exo.4
                               gauss_hump_pe4_u0.9.std.exo.5
                               gauss_hump_pe4_u0.9.std.exo.6
                               gauss_hump_pe4_u0.9.std.exo.7
                               gauss_hump_pe4_u0.9.std.exo.8
                               gauss_hump_pe4_u0.9.std.exo.9
                               gauss_hump_pe4_u0.9.std.exo.10
                               gauss_hump_pe4_u0.9.std.exo.11
                               gauss_hump_pe4_u0.9.std.exo.12
                               gauss_hump_pe4_u0.9.std.exo.13
                               gauss_hump_pe4_u0.9.std.exo.14
                               gauss_hump_pe4_u0.9.std.exo.15
                               gauss_hump_pe4_u0.9.std.exo.16
                               gauss_hump_pe4_u0.9.std.exo.17
                               gauss_hump_pe4_u0.9.std.exo.18
                               gauss_hump_pe4_u0.9.std.exo.19
                               gauss_hump_pe4_u0.9.std.exo.20
                               gauss_hump_pe4_u0.9.std.exo.21
                               gauss_hump_pe4_u0.9.std.exo.22
                               gauss_hump_pe4_u0.9.std.exo.23
                               gauss_hump_pe4_u0.9.std.exo.24
                               gauss_hump_pe4_u0.9.std.exo.25
                               gauss_hump_pe4_u0.9.std.exo.26
                               gauss_hump_pe4_u0.9.std.exo.27
                               gauss_hump_pe4_u0.9.std.exo.28
                               gauss_hump_pe4_u0.9.std.exo.29
                               gauss_hump_pe4_u0.9.std.exo.30
                               gauss_hump_pe4_u0.9.std.exo.31
                               gauss_hump_pe4_u0.9.std.exo.32
                               gauss_hump_pe4_u0.9.std.exo.33
                               gauss_hump_pe4_u0.9.std.exo.34
                               gauss_hump_pe4_u0.9.std.exo.35
                               gauss_hump_pe4_u0.9.std.exo.36
                               gauss_hump_pe4_u0.9.std.exo.37
                               gauss_hump_pe4_u0.9.std.exo.38
                               gauss_hump_pe4_u0.9.std.exo.39
                    ARGS -c gauss_hump.q -i unitsquare_01_3.6k.exo -v -u 0.9
                         +balancer RandCentLB +LBDebug 1
                    BIN_BASELINE gauss_hump_pe4_u0.9.std.exo.0
                                 gauss_hump_pe4_u0.9.std.exo.1
                                 gauss_hump_pe4_u0.9.std.exo.2
                                 gauss_hump_pe4_u0.9.std.exo.3
                                 gauss_hump_pe4_u0.9.std.exo.4
                                 gauss_hump_pe4_u0.9.std.exo.5
                                 gauss_hump_pe4_u0.9.std.exo.6
                                 gauss_hump_pe4_u0.9.std.exo.7
                                 gauss_hump_pe4_u0.9.std.exo.8
                                 gauss_hump_pe4_u0.9.std.exo.9
                                 gauss_hump_pe4_u0.9.std.exo.10
                                 gauss_hump_pe4_u0.9.std.exo.11
                                 gauss_hump_pe4_u0.9.std.exo.12
                                 gauss_hump_pe4_u0.9.std.exo.13
                                 gauss_hump_pe4_u0.9.std.exo.14
                                 gauss_hump_pe4_u0.9.std.exo.15
                                 gauss_hump_pe4_u0.9.std.exo.16
                                 gauss_hump_pe4_u0.9.std.exo.17
                                 gauss_hump_pe4_u0.9.std.exo.18
                                 gauss_hump_pe4_u0.9.std.exo.19
                                 gauss_hump_pe4_u0.9.std.exo.20
                                 gauss_hump_pe4_u0.9.std.exo.21
                                 gauss_hump_pe4_u0.9.std.exo.22
                                 gauss_hump_pe4_u0.9.std.exo.23
                                 gauss_hump_pe4_u0.9.std.exo.24
                                 gauss_hump_pe4_u0.9.std.exo.25
                                 gauss_hump_pe4_u0.9.std.exo.26
                                 gauss_hump_pe4_u0.9.std.exo.27
                                 gauss_hump_pe4_u0.9.std.exo.28
                                 gauss_hump_pe4_u0.9.std.exo.29
                                 gauss_hump_pe4_u0.9.std.exo.30
                                 gauss_hump_pe4_u0.9.std.exo.31
                                 gauss_hump_pe4_u0.9.std.exo.32
                                 gauss_hump_pe4_u0.9.std.exo.33
                                 gauss_hump_pe4_u0.9.std.exo.34
                                 gauss_hump_pe4_u0.9.std.exo.35
                                 gauss_hump_pe4_u0.9.std.exo.36
                                 gauss_hump_pe4_u0.9.std.exo.37
                                 gauss_hump_pe4_u0.9.std.exo.38
                                 gauss_hump_pe4_u0.9.std.exo.39
                    BIN_RESULT out.0 out.1 out.2 out.3 out.4 out.5 out.6 out.7
                               out.8 out.9 out.10 out.11 out.12 out.13 out.14
                               out.15 out.16 out.17 out.18 out.19 out.20 out.21
                               out.22 out.23 out.24 out.25 out.26 out.27 out.28
                               out.29 out.30 out.31 out.32 out.33 out.34 out.35
                               out.36 out.37 out.38 out.39
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
                    LABELS migration)
