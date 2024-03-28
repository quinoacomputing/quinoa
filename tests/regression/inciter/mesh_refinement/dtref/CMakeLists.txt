# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# See regression/inciter/mesh_refinement/README.txt for naming conventions

# Serial

#add_regression_test(amr_dtref_u_trans_gauss_hump_dg ${INCITER_EXECUTABLE}
#                    NUMPES 1
#                    INPUTFILES gauss_hump.q unitcube_01_112_ss3.exo
#                    ARGS -c gauss_hump.q -i unitcube_01_112_ss3.exo -v
#                    BIN_BASELINE gauss_hump_u_trans_pe1_u0.0.std.e-s.0.1.0
#                                 gauss_hump_u_trans_pe1_u0.0.std.e-s.1.1.0
#                                 gauss_hump_u_trans_pe1_u0.0.std.e-s.2.1.0
#                    BIN_RESULT out.e-s.0.1.0 out.e-s.1.1.0 out.e-s.2.1.0
#                    BIN_DIFF_PROG_CONF exodiff_gauss_hump_dg.cfg
#                    BIN_DIFF_PROG_ARGS -m
#                    TEXT_BASELINE gauss_hump_dg.std
#                    TEXT_RESULT diag
#                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
#                    LABELS dg amr)
#
#add_regression_test(amr_dtref_u_trans_reord_gauss_hump_dg ${INCITER_EXECUTABLE}
#                    NUMPES 1
#                    INPUTFILES gauss_hump_reord.q unitcube_01_112_ss3.exo
#                    ARGS -c gauss_hump_reord.q -i unitcube_01_112_ss3.exo -v
#                    BIN_BASELINE gauss_hump_u_trans_pe1_u0.0.std.e-s.0.1.0
#                                 gauss_hump_u_trans_pe1_u0.0.std.e-s.1.1.0
#                                 gauss_hump_u_trans_pe1_u0.0.std.e-s.2.1.0
#                    BIN_RESULT out.e-s.0.1.0 out.e-s.1.1.0 out.e-s.2.1.0
#                    BIN_DIFF_PROG_CONF exodiff_gauss_hump_dg.cfg
#                    BIN_DIFF_PROG_ARGS -m
#                    TEXT_BASELINE gauss_hump_dg.std
#                    TEXT_RESULT diag
#                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
#                    LABELS dg amr)

# Parallel, no virtualization

#add_regression_test(amr_dtref_u_trans_gauss_hump_dg ${INCITER_EXECUTABLE}
#                    NUMPES 4
#                    INPUTFILES gauss_hump.q unitcube_01_112_ss3.exo
#                    ARGS -c gauss_hump.q -i unitcube_01_112_ss3.exo -v
#                    BIN_BASELINE gauss_hump_u_trans_pe4_u0.0.std.e-s.0.4.0
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.0.4.1
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.0.4.2
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.0.4.3
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.1.4.0
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.1.4.1
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.1.4.2
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.1.4.3
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.2.4.0
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.2.4.1
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.2.4.2
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.2.4.3
#                    BIN_RESULT out.e-s.0.4.0
#                               out.e-s.0.4.1
#                               out.e-s.0.4.2
#                               out.e-s.0.4.3
#                               out.e-s.1.4.0
#                               out.e-s.1.4.1
#                               out.e-s.1.4.2
#                               out.e-s.1.4.3
#                               out.e-s.2.4.0
#                               out.e-s.2.4.1
#                               out.e-s.2.4.2
#                               out.e-s.2.4.3
#                    BIN_DIFF_PROG_CONF exodiff_gauss_hump_dg.cfg
#                    BIN_DIFF_PROG_ARGS -m
#                    TEXT_BASELINE gauss_hump_dg.std
#                    TEXT_RESULT diag
#                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
#                    LABELS dg amr)
#
#add_regression_test(amr_dtref_u_trans_reord_gauss_hump_dg ${INCITER_EXECUTABLE}
#                    NUMPES 4
#                    INPUTFILES gauss_hump_reord.q unitcube_01_112_ss3.exo
#                    ARGS -c gauss_hump_reord.q -i unitcube_01_112_ss3.exo -v
#                    BIN_BASELINE gauss_hump_u_trans_pe4_u0.0.std.e-s.0.4.0
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.0.4.1
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.0.4.2
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.0.4.3
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.1.4.0
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.1.4.1
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.1.4.2
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.1.4.3
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.2.4.0
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.2.4.1
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.2.4.2
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.2.4.3
#                    BIN_RESULT out.e-s.0.4.0
#                               out.e-s.0.4.1
#                               out.e-s.0.4.2
#                               out.e-s.0.4.3
#                               out.e-s.1.4.0
#                               out.e-s.1.4.1
#                               out.e-s.1.4.2
#                               out.e-s.1.4.3
#                               out.e-s.2.4.0
#                               out.e-s.2.4.1
#                               out.e-s.2.4.2
#                               out.e-s.2.4.3
#                    BIN_DIFF_PROG_CONF exodiff_gauss_hump_dg.cfg
#                    BIN_DIFF_PROG_ARGS -m
#                    TEXT_BASELINE gauss_hump_dg.std
#                    TEXT_RESULT diag
#                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
#                    LABELS dg amr)
#
#add_regression_test(amr_dtref_u_trans_gauss_hump_dg_migr ${INCITER_EXECUTABLE}
#                    NUMPES 4
#                    INPUTFILES gauss_hump.q unitcube_01_112_ss3.exo
#                    ARGS -c gauss_hump.q -i unitcube_01_112_ss3.exo -v
#                         +balancer RandCentLB +LBDebug 1 +cs
#                    BIN_BASELINE gauss_hump_u_trans_pe4_u0.0.std.e-s.0.4.0
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.0.4.1
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.0.4.2
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.0.4.3
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.1.4.0
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.1.4.1
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.1.4.2
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.1.4.3
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.2.4.0
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.2.4.1
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.2.4.2
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.2.4.3
#                    BIN_RESULT out.e-s.0.4.0
#                               out.e-s.0.4.1
#                               out.e-s.0.4.2
#                               out.e-s.0.4.3
#                               out.e-s.1.4.0
#                               out.e-s.1.4.1
#                               out.e-s.1.4.2
#                               out.e-s.1.4.3
#                               out.e-s.2.4.0
#                               out.e-s.2.4.1
#                               out.e-s.2.4.2
#                               out.e-s.2.4.3
#                    BIN_DIFF_PROG_CONF exodiff_gauss_hump_dg.cfg
#                    BIN_DIFF_PROG_ARGS -m
#                    TEXT_BASELINE gauss_hump_dg.std
#                    TEXT_RESULT diag
#                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
#                    LABELS dg amr migration)
#
#add_regression_test(amr_dtref_u_trans_reord_gauss_hump_dg_migr
#                    ${INCITER_EXECUTABLE}
#                    NUMPES 4
#                    INPUTFILES gauss_hump_reord.q unitcube_01_112_ss3.exo
#                    ARGS -c gauss_hump_reord.q -i unitcube_01_112_ss3.exo -v
#                         +balancer RandCentLB +LBDebug 1 +cs
#                    BIN_BASELINE gauss_hump_u_trans_pe4_u0.0.std.e-s.0.4.0
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.0.4.1
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.0.4.2
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.0.4.3
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.1.4.0
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.1.4.1
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.1.4.2
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.1.4.3
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.2.4.0
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.2.4.1
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.2.4.2
#                                 gauss_hump_u_trans_pe4_u0.0.std.e-s.2.4.3
#                    BIN_RESULT out.e-s.0.4.0
#                               out.e-s.0.4.1
#                               out.e-s.0.4.2
#                               out.e-s.0.4.3
#                               out.e-s.1.4.0
#                               out.e-s.1.4.1
#                               out.e-s.1.4.2
#                               out.e-s.1.4.3
#                               out.e-s.2.4.0
#                               out.e-s.2.4.1
#                               out.e-s.2.4.2
#                               out.e-s.2.4.3
#                    BIN_DIFF_PROG_CONF exodiff_gauss_hump_dg.cfg
#                    BIN_DIFF_PROG_ARGS -m
#                    TEXT_BASELINE gauss_hump_dg.std
#                    TEXT_RESULT diag
#                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
#                    LABELS dg amr migration)

# Parallel, virtualization

#add_regression_test(amr_dtref_u_trans_gauss_hump_dg_u0.8 ${INCITER_EXECUTABLE}
#                    NUMPES 4
#                    INPUTFILES gauss_hump.q unitcube_01_112_ss3.exo
#                    ARGS -c gauss_hump.q -i unitcube_01_112_ss3.exo -v -u 0.8
#                    BIN_BASELINE gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.0
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.1
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.2
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.3
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.4
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.5
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.6
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.7
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.8
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.9
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.10
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.11
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.12
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.13
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.14
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.15
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.16
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.17
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.0
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.1
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.2
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.3
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.4
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.5
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.6
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.7
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.8
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.9
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.10
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.11
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.12
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.13
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.14
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.15
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.16
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.17
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.0
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.1
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.2
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.3
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.4
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.5
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.6
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.7
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.8
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.9
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.10
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.11
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.12
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.13
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.14
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.15
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.16
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.17
#                    BIN_RESULT out.e-s.0.18.0
#                               out.e-s.0.18.1
#                               out.e-s.0.18.2
#                               out.e-s.0.18.3
#                               out.e-s.0.18.4
#                               out.e-s.0.18.5
#                               out.e-s.0.18.6
#                               out.e-s.0.18.7
#                               out.e-s.0.18.8
#                               out.e-s.0.18.9
#                               out.e-s.0.18.10
#                               out.e-s.0.18.11
#                               out.e-s.0.18.12
#                               out.e-s.0.18.13
#                               out.e-s.0.18.14
#                               out.e-s.0.18.15
#                               out.e-s.0.18.16
#                               out.e-s.0.18.17
#                               out.e-s.1.18.0
#                               out.e-s.1.18.1
#                               out.e-s.1.18.2
#                               out.e-s.1.18.3
#                               out.e-s.1.18.4
#                               out.e-s.1.18.5
#                               out.e-s.1.18.6
#                               out.e-s.1.18.7
#                               out.e-s.1.18.8
#                               out.e-s.1.18.9
#                               out.e-s.1.18.10
#                               out.e-s.1.18.11
#                               out.e-s.1.18.12
#                               out.e-s.1.18.13
#                               out.e-s.1.18.14
#                               out.e-s.1.18.15
#                               out.e-s.1.18.16
#                               out.e-s.1.18.17
#                               out.e-s.2.18.0
#                               out.e-s.2.18.1
#                               out.e-s.2.18.2
#                               out.e-s.2.18.3
#                               out.e-s.2.18.4
#                               out.e-s.2.18.5
#                               out.e-s.2.18.6
#                               out.e-s.2.18.7
#                               out.e-s.2.18.8
#                               out.e-s.2.18.9
#                               out.e-s.2.18.10
#                               out.e-s.2.18.11
#                               out.e-s.2.18.12
#                               out.e-s.2.18.13
#                               out.e-s.2.18.14
#                               out.e-s.2.18.15
#                               out.e-s.2.18.16
#                               out.e-s.2.18.17
#                    BIN_DIFF_PROG_CONF exodiff_gauss_hump_dg.cfg
#                    BIN_DIFF_PROG_ARGS -m
#                    TEXT_BASELINE gauss_hump_dg.std
#                    TEXT_RESULT diag
#                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
#                    LABELS dg amr)

#add_regression_test(amr_dtref_u_trans_gauss_hump_dg_u0.8_migr
#                    ${INCITER_EXECUTABLE}
#                    NUMPES 4
#                    INPUTFILES gauss_hump.q unitcube_01_112_ss3.exo
#                    ARGS -c gauss_hump.q -i unitcube_01_112_ss3.exo -v -u 0.8
#                         +balancer RandCentLB +LBDebug 1 +cs
#                    BIN_BASELINE gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.0
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.1
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.2
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.3
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.4
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.5
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.6
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.7
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.8
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.9
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.10
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.11
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.12
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.13
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.14
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.15
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.16
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.0.18.17
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.0
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.1
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.2
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.3
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.4
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.5
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.6
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.7
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.8
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.9
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.10
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.11
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.12
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.13
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.14
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.15
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.16
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.1.18.17
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.0
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.1
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.2
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.3
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.4
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.5
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.6
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.7
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.8
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.9
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.10
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.11
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.12
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.13
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.14
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.15
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.16
#                                 gauss_hump_u_trans_pe4_u0.8.std.e-s.2.18.17
#                    BIN_RESULT out.e-s.0.18.0
#                               out.e-s.0.18.1
#                               out.e-s.0.18.2
#                               out.e-s.0.18.3
#                               out.e-s.0.18.4
#                               out.e-s.0.18.5
#                               out.e-s.0.18.6
#                               out.e-s.0.18.7
#                               out.e-s.0.18.8
#                               out.e-s.0.18.9
#                               out.e-s.0.18.10
#                               out.e-s.0.18.11
#                               out.e-s.0.18.12
#                               out.e-s.0.18.13
#                               out.e-s.0.18.14
#                               out.e-s.0.18.15
#                               out.e-s.0.18.16
#                               out.e-s.0.18.17
#                               out.e-s.1.18.0
#                               out.e-s.1.18.1
#                               out.e-s.1.18.2
#                               out.e-s.1.18.3
#                               out.e-s.1.18.4
#                               out.e-s.1.18.5
#                               out.e-s.1.18.6
#                               out.e-s.1.18.7
#                               out.e-s.1.18.8
#                               out.e-s.1.18.9
#                               out.e-s.1.18.10
#                               out.e-s.1.18.11
#                               out.e-s.1.18.12
#                               out.e-s.1.18.13
#                               out.e-s.1.18.14
#                               out.e-s.1.18.15
#                               out.e-s.1.18.16
#                               out.e-s.1.18.17
#                               out.e-s.2.18.0
#                               out.e-s.2.18.1
#                               out.e-s.2.18.2
#                               out.e-s.2.18.3
#                               out.e-s.2.18.4
#                               out.e-s.2.18.5
#                               out.e-s.2.18.6
#                               out.e-s.2.18.7
#                               out.e-s.2.18.8
#                               out.e-s.2.18.9
#                               out.e-s.2.18.10
#                               out.e-s.2.18.11
#                               out.e-s.2.18.12
#                               out.e-s.2.18.13
#                               out.e-s.2.18.14
#                               out.e-s.2.18.15
#                               out.e-s.2.18.16
#                               out.e-s.2.18.17
#                    BIN_DIFF_PROG_CONF exodiff_gauss_hump_dg.cfg
#                    BIN_DIFF_PROG_ARGS -m
#                    TEXT_BASELINE gauss_hump_dg.std
#                    TEXT_RESULT diag
#                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
#                    LABELS dg amr migration)
