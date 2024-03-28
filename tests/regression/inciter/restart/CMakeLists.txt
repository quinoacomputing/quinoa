# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial, checkpoint on 1 PE, restart on 1 PE

add_regression_test(trans_alecg_checkpoint ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_alecg.q unitsquare_01_3.6k.exo
                    ARGS -c slot_cyl_alecg.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE slot_cyl_alecg_checkpoint.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE slot_cyl_alecg_checkpoint.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF slot_cyl_diag.ndiff.cfg
                    LABELS alecg restart)

add_regression_test(trans_alecg_restart ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_alecg_restart.q unitsquare_01_3.6k.exo
                    ARGS -c slot_cyl_alecg_restart.q -i unitsquare_01_3.6k.exo -v
                         +restart restart
                    BIN_BASELINE slot_cyl_alecg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE slot_cyl_alecg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF slot_cyl_diag.ndiff.cfg
                    CHECKPOINT trans_alecg_checkpoint
                    LABELS alecg restart)

# Parallel, checkpoint on 4 PEs, restart on 4 PEs

add_regression_test(trans_alecg_checkpoint ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl_alecg.q unitsquare_01_3.6k.exo
                    ARGS -c slot_cyl_alecg.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE slot_cyl_alecg_checkpoint.std.exo.0
                                 slot_cyl_alecg_checkpoint.std.exo.1
                                 slot_cyl_alecg_checkpoint.std.exo.2
                                 slot_cyl_alecg_checkpoint.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE slot_cyl_alecg_checkpoint.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF slot_cyl_alecg_diag.ndiff.cfg
                    LABELS alecg restart)

add_regression_test(trans_alecg_restart ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl_alecg_restart.q unitsquare_01_3.6k.exo
                    ARGS -c slot_cyl_alecg_restart.q -i unitsquare_01_3.6k.exo -v
                         +restart restart
                    BIN_BASELINE slot_cyl_alecg.std.exo.0
                                 slot_cyl_alecg.std.exo.1
                                 slot_cyl_alecg.std.exo.2
                                 slot_cyl_alecg.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE slot_cyl_alecg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF slot_cyl_alecg_diag.ndiff.cfg
                    CHECKPOINT trans_alecg_checkpoint
                    LABELS alecg restart)
