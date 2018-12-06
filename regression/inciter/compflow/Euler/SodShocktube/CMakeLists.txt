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
