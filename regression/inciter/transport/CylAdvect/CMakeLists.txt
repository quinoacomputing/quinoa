# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Virtualization = 0

add_regression_test(cyl_advect_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES cyl_advect_dg.q unitsquare_01_3.6k.exo diag_dg.std
                               exodiff.cfg cyl_advect_dg.std.exo
                    ARGS -c cyl_advect_dg.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE cyl_advect_dg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF cyl_advect_diag.ndiff.cfg)

add_regression_test(cyl_advect_dgp1 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES cyl_advect_dgp1.q unitsquare_01_3.6k.exo diag_dgp1.std
                               exodiff.cfg cyl_advect_dgp1.std.exo
                    ARGS -c cyl_advect_dgp1.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE cyl_advect_dgp1.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dgp1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF cyl_advect_diag.ndiff.cfg)
