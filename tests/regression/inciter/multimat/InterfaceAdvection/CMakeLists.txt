# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(multimat_interface_advection_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES interface_advection_dg.lua unitsquare_01_3.6k.exo
                    ARGS -c interface_advection_dg.lua -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE interface_advection_dg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF interface_advection_diag.ndiff.cfg
                    LABELS dg multimat)

# Parallel

add_regression_test(multimat_interface_advection_dg ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES interface_advection_dg.lua unitsquare_01_3.6k.exo
                    ARGS -c interface_advection_dg.lua -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE interface_advection_dg_pe4_u0.0.std.exo.0
                                 interface_advection_dg_pe4_u0.0.std.exo.1
                                 interface_advection_dg_pe4_u0.0.std.exo.2
                                 interface_advection_dg_pe4_u0.0.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    TEXT_BASELINE diag_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF interface_advection_diag.ndiff.cfg
                    LABELS dg multimat)

add_regression_test(multimat_1dinterface_advection_dg ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES interface_1dadvection_dg.lua rectangle_01_int50_struc896.exo
                    ARGS -c interface_1dadvection_dg.lua -i rectangle_01_int50_struc896.exo -v
                    BIN_BASELINE interface_1dadvection_dg_pe4_u0.0.std.exo.0
                                 interface_1dadvection_dg_pe4_u0.0.std.exo.1
                                 interface_1dadvection_dg_pe4_u0.0.std.exo.2
                                 interface_1dadvection_dg_pe4_u0.0.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff_1ddg.cfg
                    TEXT_BASELINE diag_1ddg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF interface_advection_diag.ndiff.cfg
                    LABELS dg multimat)
