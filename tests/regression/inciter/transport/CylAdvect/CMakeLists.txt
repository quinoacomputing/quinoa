# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Virtualization = 0

add_regression_test(cyl_advect_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES cyl_advect_dg.lua unitsquare_01_3.6k.exo
                    ARGS -c cyl_advect_dg.lua -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE cyl_advect_dg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF cyl_advect_diag.ndiff.cfg
                    LABELS dg)

add_regression_test(cyl_advect_dgp1 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES cyl_advect_dgp1.lua unitsquare_01_3.6k.exo
                    ARGS -c cyl_advect_dgp1.lua -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE cyl_advect_dgp1.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dgp1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF cyl_advect_diag.ndiff.cfg
                    LABELS dg)

add_regression_test(cyl_advect_dgp1 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES cyl_advect_dgp1.lua unitsquare_01_3.6k.exo
                    ARGS -c cyl_advect_dgp1.lua -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE cyl_advect_dgp1_pe4_u0.0.std.exo.0
                                 cyl_advect_dgp1_pe4_u0.0.std.exo.1
                                 cyl_advect_dgp1_pe4_u0.0.std.exo.2
                                 cyl_advect_dgp1_pe4_u0.0.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_dgp1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF cyl_advect_diag.ndiff.cfg
                    LABELS dg)

# Virtualization = 0.90

add_regression_test(cyl_advect_dgp1_u0.90 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES cyl_advect_dgp1.lua unitsquare_01_3.6k.exo
                    ARGS -c cyl_advect_dgp1.lua -i unitsquare_01_3.6k.exo -v -u 0.90
                    BIN_BASELINE cyl_advect_dgp1_pe4_u0.90.std.exo.0
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.1
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.2
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.3
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.4
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.5
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.6
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.7
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.8
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.9
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.10
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.11
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.12
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.13
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.14
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.15
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.16
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.17
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.18
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.19
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.20
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.21
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.22
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.23
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.24
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.25
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.26
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.27
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.28
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.29
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.30
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.31
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.32
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.33
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.34
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.35
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.36
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.37
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.38
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.39
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
                    TEXT_BASELINE diag_dgp1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF cyl_advect_diag.ndiff.cfg
                    LABELS dg)

# Migration + Virtualization(0.90)

add_regression_test(cyl_advect_dgp1_u0.90_lbfreq5_migr ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES cyl_advect_dgp1.lua unitsquare_01_3.6k.exo
                    ARGS -c cyl_advect_dgp1.lua -i unitsquare_01_3.6k.exo -v -u 0.90
                         +balancer RandCentLB +LBDebug 1 -l 5
                    BIN_BASELINE cyl_advect_dgp1_pe4_u0.90.std.exo.0
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.1
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.2
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.3
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.4
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.5
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.6
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.7
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.8
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.9
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.10
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.11
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.12
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.13
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.14
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.15
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.16
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.17
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.18
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.19
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.20
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.21
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.22
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.23
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.24
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.25
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.26
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.27
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.28
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.29
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.30
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.31
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.32
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.33
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.34
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.35
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.36
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.37
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.38
                                 cyl_advect_dgp1_pe4_u0.90.std.exo.39
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
                    TEXT_BASELINE diag_dgp1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF cyl_advect_diag.ndiff.cfg
                    LABELS dg migration)

add_regression_test(cyl_advect_dgp1_weno_u0.90_lbfreq5_migr ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES cyl_advect_dgp1_weno.lua unitsquare_01_3.6k.exo
                    ARGS -c cyl_advect_dgp1_weno.lua -i unitsquare_01_3.6k.exo -v -u 0.90
                         +balancer RandCentLB +LBDebug 1 -l 5
                    BIN_BASELINE cyl_advect_dgp1_weno_pe4_u0.90.std.exo.0
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.1
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.2
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.3
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.4
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.5
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.6
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.7
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.8
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.9
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.10
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.11
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.12
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.13
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.14
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.15
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.16
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.17
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.18
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.19
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.20
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.21
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.22
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.23
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.24
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.25
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.26
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.27
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.28
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.29
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.30
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.31
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.32
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.33
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.34
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.35
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.36
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.37
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.38
                                 cyl_advect_dgp1_weno_pe4_u0.90.std.exo.39
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
                    TEXT_BASELINE diag_dgp1_weno.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF cyl_advect_diag.ndiff.cfg
                    LABELS dg migration)
