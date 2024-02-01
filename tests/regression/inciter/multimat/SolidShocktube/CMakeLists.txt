# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Parallel
add_regression_test(multimat_solid_shocktube_p0p1 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES solid_shocktube_p0p1.q rectangle_01_int50_struc896.exo
                    ARGS -c solid_shocktube_p0p1.q -i rectangle_01_int50_struc896.exo -v
                    BIN_BASELINE solid_shocktube_p0p1.std.exo.0
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_p0p1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF ndiff.cfg
                    LABELS dg multimat solids)
