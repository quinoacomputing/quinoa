# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(multimat_sod_shocktube_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod_shocktube_dg.q rectangle_01_1.5k.exo
                    ARGS -c sod_shocktube_dg.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod_shocktube_dg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    TEXT_BASELINE diag_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg
                    LABELS dg multimat)

add_regression_test(multimat_sod_shocktube_p0p1 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod_shocktube_p0p1.q rectangle_01_1.5k.exo
                    ARGS -c sod_shocktube_p0p1.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod_shocktube_p0p1.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    TEXT_BASELINE diag_p0p1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg
                    LABELS dg multimat)

add_regression_test(multimat_sod_shocktube_dgp1 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod_shocktube_dgp1.q rectangle_01_1.5k.exo
                    ARGS -c sod_shocktube_dgp1.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod_shocktube_dgp1.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dgp1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg
                    LABELS dg)

add_regression_test(multimat_sod_shocktube_dgp1_limsolprojection
                    ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod_shocktube_dgp1_lsp.q rectangle_01_1.5k.exo
                    ARGS -c sod_shocktube_dgp1_lsp.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod_shocktube_dgp1_lsp.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dgp1_lsp.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg
                    LABELS dg)

add_regression_test(multimat_sod_shocktube_fv_hll ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod_shocktube_fv_hll.q rectangle_01_1.5k.exo
                    ARGS -c sod_shocktube_fv_hll.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod_shocktube_fv_hll.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_fv.cfg
                    TEXT_BASELINE diag_fv_hll.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg
                    LABELS fv multimat)

# Parallel

add_regression_test(multimat_sod_shocktube_dg ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES sod_shocktube_dg.q rectangle_01_1.5k.exo
                    ARGS -c sod_shocktube_dg.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod_shocktube_dg_pe4_u0.0.std.exo.0
                                 sod_shocktube_dg_pe4_u0.0.std.exo.1
                                 sod_shocktube_dg_pe4_u0.0.std.exo.2
                                 sod_shocktube_dg_pe4_u0.0.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    TEXT_BASELINE diag_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg
                    LABELS dg multimat)

add_regression_test(multimat_sod_shocktube_p0p1 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES sod_shocktube_p0p1.q rectangle_01_1.5k.exo
                    ARGS -c sod_shocktube_p0p1.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod_shocktube_p0p1_pe4_u0.0.std.exo.0
                                 sod_shocktube_p0p1_pe4_u0.0.std.exo.1
                                 sod_shocktube_p0p1_pe4_u0.0.std.exo.2
                                 sod_shocktube_p0p1_pe4_u0.0.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    TEXT_BASELINE diag_p0p1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg
                    LABELS dg multimat)

add_regression_test(multimat_sod_shocktube_dgp1 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES sod_shocktube_dgp1.q rectangle_01_1.5k.exo
                    ARGS -c sod_shocktube_dgp1.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE sod_shocktube_dgp1_pe4_u0.0.std.exo.0
                                 sod_shocktube_dgp1_pe4_u0.0.std.exo.1
                                 sod_shocktube_dgp1_pe4_u0.0.std.exo.2
                                 sod_shocktube_dgp1_pe4_u0.0.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    TEXT_BASELINE diag_dgp1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg
                    LABELS dg)

# Parallel + Virtualization

add_regression_test(multimat_sod_shocktube_p0p1_u0.9 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES sod_shocktube_p0p1.q rectangle_01_1.5k.exo
                    ARGS -c sod_shocktube_p0p1.q -i rectangle_01_1.5k.exo -u 0.9 -v
                    BIN_BASELINE sod_shocktube_p0p1_pe4_u0.9.std.exo.0
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.1
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.2
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.3
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.4
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.5
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.6
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.7
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.8
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.9
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.10
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.11
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.12
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.13
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.14
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.15
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.16
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.17
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.18
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.19
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.20
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.21
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.22
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.23
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.24
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.25
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.26
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.27
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.28
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.29
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.30
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.31
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.32
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.33
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.34
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.35
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.36
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.37
                                 sod_shocktube_p0p1_pe4_u0.9.std.exo.38
                    BIN_RESULT out.e-s.0.39.0
                               out.e-s.0.39.1
                               out.e-s.0.39.2
                               out.e-s.0.39.3
                               out.e-s.0.39.4
                               out.e-s.0.39.5
                               out.e-s.0.39.6
                               out.e-s.0.39.7
                               out.e-s.0.39.8
                               out.e-s.0.39.9
                               out.e-s.0.39.10
                               out.e-s.0.39.11
                               out.e-s.0.39.12
                               out.e-s.0.39.13
                               out.e-s.0.39.14
                               out.e-s.0.39.15
                               out.e-s.0.39.16
                               out.e-s.0.39.17
                               out.e-s.0.39.18
                               out.e-s.0.39.19
                               out.e-s.0.39.20
                               out.e-s.0.39.21
                               out.e-s.0.39.22
                               out.e-s.0.39.23
                               out.e-s.0.39.24
                               out.e-s.0.39.25
                               out.e-s.0.39.26
                               out.e-s.0.39.27
                               out.e-s.0.39.28
                               out.e-s.0.39.29
                               out.e-s.0.39.30
                               out.e-s.0.39.31
                               out.e-s.0.39.32
                               out.e-s.0.39.33
                               out.e-s.0.39.34
                               out.e-s.0.39.35
                               out.e-s.0.39.36
                               out.e-s.0.39.37
                               out.e-s.0.39.38
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    TEXT_BASELINE diag_p0p1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg
                    LABELS dg multimat)
