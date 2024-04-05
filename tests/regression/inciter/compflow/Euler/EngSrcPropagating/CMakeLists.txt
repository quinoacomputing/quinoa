# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Parallel

add_regression_test(compflow_engsrc_propagating_alecg
                    ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES engsrc_propagating_alecg.lua tubetest_12k.exo
                    ARGS -c engsrc_propagating_alecg.lua -i tubetest_12k.exo -v
                    BIN_BASELINE engsrc_propagating_alecg_pe4.std.exo.0
                                 engsrc_propagating_alecg_pe4.std.exo.1
                                 engsrc_propagating_alecg_pe4.std.exo.2
                                 engsrc_propagating_alecg_pe4.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_alecg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_shocktube_diag.ndiff.cfg
                    LABELS alecg)
