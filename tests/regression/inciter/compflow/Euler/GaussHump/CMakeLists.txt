# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(compflow_euler_gauss_hump_pdg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES gauss_hump_pdg.q unitsquare_01_3.6k.exo
                    ARGS -c gauss_hump_pdg.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE gauss_hump_pdg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_pdg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)

# Parallel + no virtualization

add_regression_test(compflow_euler_gauss_hump_pdg ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gauss_hump_pdg.q unitsquare_01_3.6k.exo
                    ARGS -c gauss_hump_pdg.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE gauss_hump_pdg_pe4.std.exo.0
                                 gauss_hump_pdg_pe4.std.exo.1
                                 gauss_hump_pdg_pe4.std.exo.2
                                 gauss_hump_pdg_pe4.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_pdg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)

# Parallel + virtualization

add_regression_test(compflow_euler_gauss_hump_pdg_u0.8 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gauss_hump_pdg.q unitsquare_01_3.6k.exo
                    ARGS -c gauss_hump_pdg.q -i unitsquare_01_3.6k.exo -v -u 0.8
                    BIN_BASELINE gauss_hump_pdg_pe4_u0.8.std.exo.0
                                 gauss_hump_pdg_pe4_u0.8.std.exo.1
                                 gauss_hump_pdg_pe4_u0.8.std.exo.2
                                 gauss_hump_pdg_pe4_u0.8.std.exo.3
                                 gauss_hump_pdg_pe4_u0.8.std.exo.4
                                 gauss_hump_pdg_pe4_u0.8.std.exo.5
                                 gauss_hump_pdg_pe4_u0.8.std.exo.6
                                 gauss_hump_pdg_pe4_u0.8.std.exo.7
                                 gauss_hump_pdg_pe4_u0.8.std.exo.8
                                 gauss_hump_pdg_pe4_u0.8.std.exo.9
                                 gauss_hump_pdg_pe4_u0.8.std.exo.10
                                 gauss_hump_pdg_pe4_u0.8.std.exo.11
                                 gauss_hump_pdg_pe4_u0.8.std.exo.12
                                 gauss_hump_pdg_pe4_u0.8.std.exo.13
                                 gauss_hump_pdg_pe4_u0.8.std.exo.14
                                 gauss_hump_pdg_pe4_u0.8.std.exo.15
                                 gauss_hump_pdg_pe4_u0.8.std.exo.16
                                 gauss_hump_pdg_pe4_u0.8.std.exo.17
                                 gauss_hump_pdg_pe4_u0.8.std.exo.18
                                 gauss_hump_pdg_pe4_u0.8.std.exo.19
                    BIN_RESULT out.e-s.0.20.0
                               out.e-s.0.20.1
                               out.e-s.0.20.2
                               out.e-s.0.20.3
                               out.e-s.0.20.4
                               out.e-s.0.20.5
                               out.e-s.0.20.6
                               out.e-s.0.20.7
                               out.e-s.0.20.8
                               out.e-s.0.20.9
                               out.e-s.0.20.10
                               out.e-s.0.20.11
                               out.e-s.0.20.12
                               out.e-s.0.20.13
                               out.e-s.0.20.14
                               out.e-s.0.20.15
                               out.e-s.0.20.16
                               out.e-s.0.20.17
                               out.e-s.0.20.18
                               out.e-s.0.20.19
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_pdg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)


# Parallel + virtualization + migration

add_regression_test(compflow_euler_gauss_hump_pdg_u0.8_migr ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gauss_hump_pdg.q unitsquare_01_3.6k.exo
                    ARGS -c gauss_hump_pdg.q -i unitsquare_01_3.6k.exo -v -u 0.8
                         +balancer GreedyRefineLB +LBDebug 1
                    BIN_BASELINE gauss_hump_pdg_pe4_u0.8.std.exo.0
                                 gauss_hump_pdg_pe4_u0.8.std.exo.1
                                 gauss_hump_pdg_pe4_u0.8.std.exo.2
                                 gauss_hump_pdg_pe4_u0.8.std.exo.3
                                 gauss_hump_pdg_pe4_u0.8.std.exo.4
                                 gauss_hump_pdg_pe4_u0.8.std.exo.5
                                 gauss_hump_pdg_pe4_u0.8.std.exo.6
                                 gauss_hump_pdg_pe4_u0.8.std.exo.7
                                 gauss_hump_pdg_pe4_u0.8.std.exo.8
                                 gauss_hump_pdg_pe4_u0.8.std.exo.9
                                 gauss_hump_pdg_pe4_u0.8.std.exo.10
                                 gauss_hump_pdg_pe4_u0.8.std.exo.11
                                 gauss_hump_pdg_pe4_u0.8.std.exo.12
                                 gauss_hump_pdg_pe4_u0.8.std.exo.13
                                 gauss_hump_pdg_pe4_u0.8.std.exo.14
                                 gauss_hump_pdg_pe4_u0.8.std.exo.15
                                 gauss_hump_pdg_pe4_u0.8.std.exo.16
                                 gauss_hump_pdg_pe4_u0.8.std.exo.17
                                 gauss_hump_pdg_pe4_u0.8.std.exo.18
                                 gauss_hump_pdg_pe4_u0.8.std.exo.19
                    BIN_RESULT out.e-s.0.20.0
                               out.e-s.0.20.1
                               out.e-s.0.20.2
                               out.e-s.0.20.3
                               out.e-s.0.20.4
                               out.e-s.0.20.5
                               out.e-s.0.20.6
                               out.e-s.0.20.7
                               out.e-s.0.20.8
                               out.e-s.0.20.9
                               out.e-s.0.20.10
                               out.e-s.0.20.11
                               out.e-s.0.20.12
                               out.e-s.0.20.13
                               out.e-s.0.20.14
                               out.e-s.0.20.15
                               out.e-s.0.20.16
                               out.e-s.0.20.17
                               out.e-s.0.20.18
                               out.e-s.0.20.19
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_pdg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
                    LABELS migration)
