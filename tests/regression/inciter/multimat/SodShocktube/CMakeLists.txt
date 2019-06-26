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
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg)

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
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg)
