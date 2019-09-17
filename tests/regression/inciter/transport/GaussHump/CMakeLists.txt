# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Virtualization = 0

add_regression_test(gauss_hump_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES gauss_hump.q unitsquare_01_3.6k.exo
                    ARGS -c gauss_hump.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE gauss_hump.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)

add_regression_test(gauss_hump_p0p1 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES gauss_hump_p0p1.q unitsquare_01_3.6k.exo
                    ARGS -c gauss_hump_p0p1.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE gauss_hump_p0p1.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_p0p1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)

add_regression_test(gauss_hump_dgp1 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES gauss_hump_dgp1.q unitsquare_01_3.6k.exo
                    ARGS -c gauss_hump_dgp1.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE gauss_hump_dgp1.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dgp1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)

add_regression_test(gauss_hump_dgp2 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES gauss_hump_dgp2.q unitsquare_01_3.6k.exo
                    ARGS -c gauss_hump_dgp2.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE gauss_hump_dgp2.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dgp2.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)

add_regression_test(gauss_hump_dg ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gauss_hump.q unitsquare_01_3.6k.exo
                    ARGS -c gauss_hump.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE gauss_hump_pe4_u0.0.std.exo.0
                                 gauss_hump_pe4_u0.0.std.exo.1
                                 gauss_hump_pe4_u0.0.std.exo.2
                                 gauss_hump_pe4_u0.0.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)

add_regression_test(gauss_hump_p0p1 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gauss_hump_p0p1.q unitsquare_01_3.6k.exo
                    ARGS -c gauss_hump_p0p1.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE gauss_hump_p0p1_pe4_u0.0.std.exo.0
                                 gauss_hump_p0p1_pe4_u0.0.std.exo.1
                                 gauss_hump_p0p1_pe4_u0.0.std.exo.2
                                 gauss_hump_p0p1_pe4_u0.0.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_p0p1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)

add_regression_test(gauss_hump_dgp1 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gauss_hump_dgp1.q unitsquare_01_3.6k.exo
                    ARGS -c gauss_hump_dgp1.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE gauss_hump_dgp1_pe4_u0.0.std.exo.0
                                 gauss_hump_dgp1_pe4_u0.0.std.exo.1
                                 gauss_hump_dgp1_pe4_u0.0.std.exo.2
                                 gauss_hump_dgp1_pe4_u0.0.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_dgp1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)

add_regression_test(gauss_hump_dgp2 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gauss_hump_dgp2.q unitsquare_01_3.6k.exo
                    ARGS -c gauss_hump_dgp2.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE gauss_hump_dgp2_pe4.std.exo.0
                                 gauss_hump_dgp2_pe4.std.exo.1
                                 gauss_hump_dgp2_pe4.std.exo.2
                                 gauss_hump_dgp2_pe4.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_dgp2.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)

# Parallel + virtualization

add_regression_test(gauss_hump_dg_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gauss_hump.q unitsquare_01_3.6k.exo
                    ARGS -c gauss_hump.q -i unitsquare_01_3.6k.exo -v -u 0.5
                    BIN_BASELINE gauss_hump_pe4_u0.5.std.exo.0
                                 gauss_hump_pe4_u0.5.std.exo.1
                                 gauss_hump_pe4_u0.5.std.exo.2
                                 gauss_hump_pe4_u0.5.std.exo.3
                                 gauss_hump_pe4_u0.5.std.exo.4
                                 gauss_hump_pe4_u0.5.std.exo.5
                                 gauss_hump_pe4_u0.5.std.exo.6
                                 gauss_hump_pe4_u0.5.std.exo.7
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
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)

add_regression_test(gauss_hump_dg_u0.9 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gauss_hump_cube.q unitcube_1k.exo diag_cube.std
                    ARGS -c gauss_hump_cube.q -i unitcube_1k.exo -v -u 0.9
                    TEXT_BASELINE diag_cube.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)

add_regression_test(gauss_hump_p0p1_u0.9 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gauss_hump_p0p1.q unitsquare_01_3.6k.exo
                    ARGS -c gauss_hump_p0p1.q -i unitsquare_01_3.6k.exo -v -u 0.9
                    BIN_BASELINE gauss_hump_p0p1_pe4_u0.9.std.exo.0
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.1
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.2
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.3
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.4
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.5
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.6
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.7
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.8
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.9
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.10
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.11
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.12
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.13
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.14
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.15
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.16
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.17
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.18
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.19
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.20
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.21
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.22
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.23
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.24
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.25
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.26
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.27
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.28
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.29
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.30
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.31
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.32
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.33
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.34
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.35
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.36
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.37
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.38
                                 gauss_hump_p0p1_pe4_u0.9.std.exo.39
                    BIN_RESULT out.e-s.0.40.0
                               out.e-s.0.40.1
                               out.e-s.0.40.2
                               out.e-s.0.40.3
                               out.e-s.0.40.4
                               out.e-s.0.40.5
                               out.e-s.0.40.6
                               out.e-s.0.40.7
                               out.e-s.0.40.8
                               out.e-s.0.40.9
                               out.e-s.0.40.10
                               out.e-s.0.40.11
                               out.e-s.0.40.12
                               out.e-s.0.40.13
                               out.e-s.0.40.14
                               out.e-s.0.40.15
                               out.e-s.0.40.16
                               out.e-s.0.40.17
                               out.e-s.0.40.18
                               out.e-s.0.40.19
                               out.e-s.0.40.20
                               out.e-s.0.40.21
                               out.e-s.0.40.22
                               out.e-s.0.40.23
                               out.e-s.0.40.24
                               out.e-s.0.40.25
                               out.e-s.0.40.26
                               out.e-s.0.40.27
                               out.e-s.0.40.28
                               out.e-s.0.40.29
                               out.e-s.0.40.30
                               out.e-s.0.40.31
                               out.e-s.0.40.32
                               out.e-s.0.40.33
                               out.e-s.0.40.34
                               out.e-s.0.40.35
                               out.e-s.0.40.36
                               out.e-s.0.40.37
                               out.e-s.0.40.38
                               out.e-s.0.40.39
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_p0p1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg)


# Parallel + virtualization + migration

add_regression_test(gauss_hump_dg_u0.9_migr ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gauss_hump.q unitsquare_01_3.6k.exo
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
                    BIN_RESULT out.e-s.0.40.0
                               out.e-s.0.40.1
                               out.e-s.0.40.2
                               out.e-s.0.40.3
                               out.e-s.0.40.4
                               out.e-s.0.40.5
                               out.e-s.0.40.6
                               out.e-s.0.40.7
                               out.e-s.0.40.8
                               out.e-s.0.40.9
                               out.e-s.0.40.10
                               out.e-s.0.40.11
                               out.e-s.0.40.12
                               out.e-s.0.40.13
                               out.e-s.0.40.14
                               out.e-s.0.40.15
                               out.e-s.0.40.16
                               out.e-s.0.40.17
                               out.e-s.0.40.18
                               out.e-s.0.40.19
                               out.e-s.0.40.20
                               out.e-s.0.40.21
                               out.e-s.0.40.22
                               out.e-s.0.40.23
                               out.e-s.0.40.24
                               out.e-s.0.40.25
                               out.e-s.0.40.26
                               out.e-s.0.40.27
                               out.e-s.0.40.28
                               out.e-s.0.40.29
                               out.e-s.0.40.30
                               out.e-s.0.40.31
                               out.e-s.0.40.32
                               out.e-s.0.40.33
                               out.e-s.0.40.34
                               out.e-s.0.40.35
                               out.e-s.0.40.36
                               out.e-s.0.40.37
                               out.e-s.0.40.38
                               out.e-s.0.40.39
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
                    LABELS migration)
