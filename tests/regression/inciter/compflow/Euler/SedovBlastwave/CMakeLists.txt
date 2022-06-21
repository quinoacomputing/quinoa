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
                    TEXT_DIFF_PROG_CONF sedov_blastwave_diag.ndiff.cfg
                    LABELS dg)

#add_regression_test(compflow_euler_sedovblastwave_dgp2 ${INCITER_EXECUTABLE}
#                    NUMPES 1
#                    INPUTFILES sedov_blastwave_dgp2.q unitsquare_01_3.6k.exo
#                    ARGS -c sedov_blastwave_dgp2.q -i unitsquare_01_3.6k.exo -v
#                    BIN_BASELINE sedov_blastwave_dgp2.std.exo
#                    BIN_RESULT out.e-s.0.1.0
#                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
#                    BIN_DIFF_PROG_ARGS -m
#                    TEXT_BASELINE diag_dgp2.std
#                    TEXT_RESULT diag
#                    TEXT_DIFF_PROG_CONF sedov_blastwave_diag.ndiff.cfg
#                    LABELS dg)

add_regression_test(compflow_euler_sedovblastwave_pdg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sedov_blastwave_pdg.q unitsquare_01_3.6k.exo
                    ARGS -c sedov_blastwave_pdg.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE sedov_blastwave_pdg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_pdg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sedov_blastwave_diag.ndiff.cfg
                    LABELS dg)

# Parallel

add_regression_test(compflow_euler_sedovblastwave_dgp1 ${INCITER_EXECUTABLE}
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
                    TEXT_DIFF_PROG_CONF sedov_blastwave_diag.ndiff.cfg
                    LABELS dg)

#add_regression_test(compflow_euler_sedovblastwave_dgp2 ${INCITER_EXECUTABLE}
#                    NUMPES 4
#                    INPUTFILES sedov_blastwave_dgp2.q unitsquare_01_3.6k.exo
#                    ARGS -c sedov_blastwave_dgp2.q -i unitsquare_01_3.6k.exo -v
#                    BIN_BASELINE sedov_blastwave_dgp2_pe4.std.exo.0
#                                 sedov_blastwave_dgp2_pe4.std.exo.1
#                                 sedov_blastwave_dgp2_pe4.std.exo.2
#                                 sedov_blastwave_dgp2_pe4.std.exo.3
#                    BIN_RESULT out.e-s.0.4.0
#                               out.e-s.0.4.1
#                               out.e-s.0.4.2
#                               out.e-s.0.4.3
#                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
#                    BIN_DIFF_PROG_ARGS -m
#                    TEXT_BASELINE diag_dgp2.std
#                    TEXT_RESULT diag
#                    TEXT_DIFF_PROG_CONF sedov_blastwave_diag.ndiff.cfg
#                    LABELS dg)

add_regression_test(compflow_euler_sedovblastwave_pdg ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES sedov_blastwave_pdg.q unitsquare_01_3.6k.exo
                    ARGS -c sedov_blastwave_pdg.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE sedov_blastwave_pdg_pe4_u0.0.std.exo.0
                                 sedov_blastwave_pdg_pe4_u0.0.std.exo.1
                                 sedov_blastwave_pdg_pe4_u0.0.std.exo.2
                                 sedov_blastwave_pdg_pe4_u0.0.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_pdg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sedov_blastwave_diag.ndiff.cfg
                    LABELS dg)

# Parallel + Virtualization + Migration

add_regression_test(compflow_euler_sedovblastwave_pdg_u0.9_migr ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES sedov_blastwave_pdg.q unitsquare_01_3.6k.exo
                    ARGS -c sedov_blastwave_pdg.q -i unitsquare_01_3.6k.exo -v -u 0.9
                         +balancer RandCentLB +LBDebug 1
                    BIN_BASELINE sedov_blastwave_pdg_pe4_u0.9.std.exo.0
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.1
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.2
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.3
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.4
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.5
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.6
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.7
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.8
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.9
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.10
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.11
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.12
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.13
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.14
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.15
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.16
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.17
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.18
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.19
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.20
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.21
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.22
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.23
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.24
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.25
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.26
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.27
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.28
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.29
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.30
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.31
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.32
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.33
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.34
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.35
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.36
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.37
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.38
                                 sedov_blastwave_pdg_pe4_u0.9.std.exo.39
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
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_pdg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sedov_blastwave_diag.ndiff.cfg
                    LABELS dg migration)
