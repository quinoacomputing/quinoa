# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(trans_diagcg_checkpoint ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl.q unitsquare_01_3.6k.exo
                    ARGS -c slot_cyl.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE slot_cyl_checkpoint.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE slot_cyl_checkpoint.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF slot_cyl_diag.ndiff.cfg
                    LABELS restart)

add_regression_test(trans_diagcg_restart ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_restart.q unitsquare_01_3.6k.exo
                    ARGS -c slot_cyl_restart.q -i unitsquare_01_3.6k.exo -v
                         +restart restart
                    BIN_BASELINE slot_cyl.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE slot_cyl.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF slot_cyl_diag.ndiff.cfg
                    CHECKPOINT trans_diagcg_checkpoint
                    LABELS restart)
