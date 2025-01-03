# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(multispec_tpg_shock_1D_p0p1 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES tpg_shock_1D_dirBCs_regr.lua rectangle_01_1.5k.exo
                    ARGS -c tpg_shock_1D_dirBCs_regr.lua -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE tpg_shock_1D_dirBCs_regr.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    TEXT_BASELINE diag_p0p1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF tpg_shock_1D.ndiff.cfg
                    LABELS dg multispec)

# Parallel

add_regression_test(multispec_tpg_shock_1D_p0p1 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES tpg_shock_1D_dirBCs_regr.lua rectangle_01_1.5k.exo
                    ARGS -c tpg_shock_1D_dirBCs_regr.lua -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE tpg_shock_1D_dirBCs_regr_pe4_u0.0.std.exo.0
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.0.std.exo.1
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.0.std.exo.2
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.0.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    TEXT_BASELINE diag_p0p1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF tpg_shock_1D.ndiff.cfg
                    LABELS dg multispec)

# Parallel + Virtualization

add_regression_test(multispec_tpg_shock_1D_p0p1_u0.9 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES tpg_shock_1D_dirBCs_regr.lua rectangle_01_1.5k.exo
                    ARGS -c tpg_shock_1D_dirBCs_regr.lua -i rectangle_01_1.5k.exo -u 0.9 -v
                    BIN_BASELINE tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.0
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.1
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.2
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.3
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.4
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.5
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.6
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.7
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.8
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.9
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.10
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.11
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.12
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.13
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.14
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.15
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.16
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.17
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.18
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.19
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.20
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.21
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.22
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.23
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.24
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.25
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.26
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.27
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.28
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.29
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.30
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.31
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.32
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.33
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.34
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.35
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.36
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.37
                                 tpg_shock_1D_dirBCs_regr_pe4_u0.9.std.exo.38
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
                    TEXT_DIFF_PROG_CONF tpg_shock_1D.ndiff.cfg
                    LABELS dg multispec)
