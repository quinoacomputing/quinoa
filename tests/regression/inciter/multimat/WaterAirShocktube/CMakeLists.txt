# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(multimat_waterair_shocktube_p0p1 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES waterair_shocktube_p0p1.q rectangle_01_1.5k.exo
                    ARGS -c waterair_shocktube_p0p1.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE waterair_shocktube_p0p1.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_p0p1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF waterair_shocktube_diag.ndiff.cfg
                    LABELS dg)

add_regression_test(multimat_waterair_shocktube_dgp1 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES waterair_shocktube_dgp1.q rectangle_01_1.5k.exo
                    ARGS -c waterair_shocktube_dgp1.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE waterair_shocktube_dgp1.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dgp1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF waterair_shocktube_diag.ndiff.cfg
                    LABELS dg)

# Parallel

add_regression_test(multimat_waterair_shocktube_dgp1 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES waterair_shocktube_dgp1.q rectangle_01_1.5k.exo
                    ARGS -c waterair_shocktube_dgp1.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE waterair_shocktube_dgp1_pe4.std.exo.0
                                 waterair_shocktube_dgp1_pe4.std.exo.1
                                 waterair_shocktube_dgp1_pe4.std.exo.2
                                 waterair_shocktube_dgp1_pe4.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dgp1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF waterair_shocktube_diag.ndiff.cfg
                    LABELS dg)
