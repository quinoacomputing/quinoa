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

add_regression_test(overset_checkpoint ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES sphereflow_overset.lua freestream_BGmesh_15k.exo sphere_OSmesh_12k.exo
                    ARGS -c sphereflow_overset.lua -v
                    BIN_BASELINE background_checkpoint_pe4.std.exo.0
                                 background_checkpoint_pe4.std.exo.1
                                 background_checkpoint_pe4.std.exo.2
                                 background_checkpoint_pe4.std.exo.3
                                 overset_checkpoint_pe4.std.exo.0
                                 overset_checkpoint_pe4.std.exo.1
                                 overset_checkpoint_pe4.std.exo.2
                                 overset_checkpoint_pe4.std.exo.3
                    BIN_RESULT out.0.e-s.0.4.0
                               out.0.e-s.0.4.1
                               out.0.e-s.0.4.2
                               out.0.e-s.0.4.3
                               out.1.e-s.0.4.0
                               out.1.e-s.0.4.1
                               out.1.e-s.0.4.2
                               out.1.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff_overset.cfg
                    BIN_DIFF_PROG_ARGS -m
                    LABELS oversetfe restart)

add_regression_test(overset_restart ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES sphereflow_overset_restart.lua freestream_BGmesh_15k.exo sphere_OSmesh_12k.exo
                    ARGS -c sphereflow_overset_restart.lua -v +restart restart
                    BIN_BASELINE background_pe4.std.exo.0
                                 background_pe4.std.exo.1
                                 background_pe4.std.exo.2
                                 background_pe4.std.exo.3
                                 overset_pe4.std.exo.0
                                 overset_pe4.std.exo.1
                                 overset_pe4.std.exo.2
                                 overset_pe4.std.exo.3
                    BIN_RESULT out.0.e-s.0.4.0
                               out.0.e-s.0.4.1
                               out.0.e-s.0.4.2
                               out.0.e-s.0.4.3
                               out.1.e-s.0.4.0
                               out.1.e-s.0.4.1
                               out.1.e-s.0.4.2
                               out.1.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff_overset.cfg
                    BIN_DIFF_PROG_ARGS -m
                    CHECKPOINT overset_checkpoint
                    LABELS oversetfe restart)
