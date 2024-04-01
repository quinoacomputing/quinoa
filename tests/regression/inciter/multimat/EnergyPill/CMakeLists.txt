# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Parallel

add_regression_test(multimat_energy_pill_fv
                    ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES energy_pill_fv.lua tubetest_2blk_924.exo
                    ARGS -c energy_pill_fv.lua -i tubetest_2blk_924.exo -v
                    BIN_BASELINE energy_pill_fv_pe4.std.exo.0
                                 energy_pill_fv_pe4.std.exo.1
                                 energy_pill_fv_pe4.std.exo.2
                                 energy_pill_fv_pe4.std.exo.3
                                 energy_pill_fv_pe4-surf.1.std.exo.0
                                 energy_pill_fv_pe4-surf.1.std.exo.1
                                 energy_pill_fv_pe4-surf.1.std.exo.2
                                 energy_pill_fv_pe4-surf.1.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                               out-surf.1.e-s.0.4.0
                               out-surf.1.e-s.0.4.1
                               out-surf.1.e-s.0.4.2
                               out-surf.1.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff_fv.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_fv.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF energy_pill.ndiff.cfg
                    LABELS fv multimat)
