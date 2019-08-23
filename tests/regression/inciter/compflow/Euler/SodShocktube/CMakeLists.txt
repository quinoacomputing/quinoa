# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(compflow_euler_sodshocktube_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod_shocktube_dg.q rectangle_01_1.5k.exo
                    ARGS -c sod_shocktube_dg.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod_shocktube_dg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg)

add_regression_test(compflow_euler_sodshocktube_p0p1 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod_shocktube_p0p1.q rectangle_01_1.5k.exo
                    ARGS -c sod_shocktube_p0p1.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod_shocktube_p0p1.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_p0p1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg)

add_regression_test(compflow_euler_rotated_sodshocktube_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES rotated_sod_shocktube_dg.q
                               rectangle_01_1.5k_rotated.exo
                               rotated_sod_shocktube_dg.std.exo
                               diag_rotated_dg.std exodiff_dg.cfg
                    ARGS -c rotated_sod_shocktube_dg.q
                         -i rectangle_01_1.5k_rotated.exo -v
                    BIN_BASELINE rotated_sod_shocktube_dg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_rotated_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg)

add_regression_test(compflow_euler_sodshocktube_pdg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod_shocktube_pdg.q rectangle_01_1.5k.exo
                    ARGS -c sod_shocktube_pdg.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod_shocktube_pdg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_pdg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg)

# Parallel

add_regression_test(compflow_euler_sodshocktube_pdg ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES sod_shocktube_pdg.q rectangle_01_1.5k.exo
                    ARGS -c sod_shocktube_pdg.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod_shocktube_pdg_pe4.std.exo.0
                                 sod_shocktube_pdg_pe4.std.exo.1
                                 sod_shocktube_pdg_pe4.std.exo.2
                                 sod_shocktube_pdg_pe4.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_pdg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg)

# Parallel + virtualization

add_regression_test(compflow_euler_sodshocktube_pdg_u0.8 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES sod_shocktube_pdg.q rectangle_01_1.5k.exo
                    ARGS -c sod_shocktube_pdg.q -i rectangle_01_1.5k.exo -v -u 0.8
                    BIN_BASELINE sod_shocktube_pdg_pe4_u0.8.std.exo.0
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.1
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.2
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.3
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.4
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.5
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.6
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.7
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.8
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.9
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.10
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.11
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.12
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.13
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.14
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.15
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.16
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.17
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.18
                    BIN_RESULT out.e-s.0.19.0
                               out.e-s.0.19.1
                               out.e-s.0.19.2
                               out.e-s.0.19.3
                               out.e-s.0.19.4
                               out.e-s.0.19.5
                               out.e-s.0.19.6
                               out.e-s.0.19.7
                               out.e-s.0.19.8
                               out.e-s.0.19.9
                               out.e-s.0.19.10
                               out.e-s.0.19.11
                               out.e-s.0.19.12
                               out.e-s.0.19.13
                               out.e-s.0.19.14
                               out.e-s.0.19.15
                               out.e-s.0.19.16
                               out.e-s.0.19.17
                               out.e-s.0.19.18
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_pdg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg)

# Parallel + virtualization + migration

add_regression_test(compflow_euler_sodshocktube_pdg_u0.8_migr ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES sod_shocktube_pdg.q rectangle_01_1.5k.exo
                    ARGS -c sod_shocktube_pdg.q -i rectangle_01_1.5k.exo -v -u 0.8
                    +balancer RandCentLB +LBDebug 1
                    BIN_BASELINE sod_shocktube_pdg_pe4_u0.8.std.exo.0
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.1
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.2
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.3
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.4
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.5
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.6
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.7
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.8
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.9
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.10
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.11
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.12
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.13
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.14
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.15
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.16
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.17
                                 sod_shocktube_pdg_pe4_u0.8.std.exo.18
                    BIN_RESULT out.e-s.0.19.0
                               out.e-s.0.19.1
                               out.e-s.0.19.2
                               out.e-s.0.19.3
                               out.e-s.0.19.4
                               out.e-s.0.19.5
                               out.e-s.0.19.6
                               out.e-s.0.19.7
                               out.e-s.0.19.8
                               out.e-s.0.19.9
                               out.e-s.0.19.10
                               out.e-s.0.19.11
                               out.e-s.0.19.12
                               out.e-s.0.19.13
                               out.e-s.0.19.14
                               out.e-s.0.19.15
                               out.e-s.0.19.16
                               out.e-s.0.19.17
                               out.e-s.0.19.18
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_pdg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg
                    LABELS migration)
