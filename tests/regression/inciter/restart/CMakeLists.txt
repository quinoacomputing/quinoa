# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial, checkpoint on 1 PE, restart on 1 PE

add_regression_test(multimat_fv_checkpoint ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES mmsod_fv.lua rectangle_01_1.5k.exo
                    ARGS -c mmsod_fv.lua -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE mmsod_fv_checkpoint.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE mmsod_fv_checkpoint.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF fv_diag.ndiff.cfg
                    LABELS fv restart)

add_regression_test(multimat_fv_restart ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES mmsod_fv_restart.lua rectangle_01_1.5k.exo
                    ARGS -c mmsod_fv_restart.lua -i rectangle_01_1.5k.exo -v
                         +restart restart
                    BIN_BASELINE mmsod_fv.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE mmsod_fv.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF fv_diag.ndiff.cfg
                    CHECKPOINT multimat_fv_checkpoint
                    LABELS fv restart)

# Parallel, checkpoint on 4 PEs, restart on 4 PEs

add_regression_test(multimat_fv_checkpoint ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES mmsod_fv.lua rectangle_01_1.5k.exo
                    ARGS -c mmsod_fv.lua -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE mmsod_fv_checkpoint.std.exo.0
                                 mmsod_fv_checkpoint.std.exo.1
                                 mmsod_fv_checkpoint.std.exo.2
                                 mmsod_fv_checkpoint.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE mmsod_fv_checkpoint.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF fv_diag.ndiff.cfg
                    LABELS fv restart)

add_regression_test(multimat_fv_restart ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES mmsod_fv_restart.lua rectangle_01_1.5k.exo
                    ARGS -c mmsod_fv_restart.lua -i rectangle_01_1.5k.exo -v
                         +restart restart
                    BIN_BASELINE mmsod_fv.std.exo.0
                                 mmsod_fv.std.exo.1
                                 mmsod_fv.std.exo.2
                                 mmsod_fv.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE mmsod_fv.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF fv_diag.ndiff.cfg
                    CHECKPOINT multimat_fv_checkpoint
                    LABELS fv restart)
