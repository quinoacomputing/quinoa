# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(compflow_euler_alecg_ale_gcl ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES gcl.q rectangle_01_1.5k.exo
                    ARGS -c gcl.q -i rectangle_01_1.5k.exo -v
                    BIN_BASELINE gcl.0.std.exo gcl.10.std.exo
                    BIN_RESULT out.e-s.0.1.0 out.e-s.10.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_gcl.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gcl_diag.ndiff.cfg
                    LABELS alecg ale)

add_regression_test(compflow_euler_alecg_ale_sod_sine ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod_shocktube_alecg_ale_sine.q
                               rectangle_01_1.5k.exo
                    ARGS -c sod_shocktube_alecg_ale_sine.q -v
                    BIN_BASELINE sod_shocktube_alecg_ale_sine.0.std.exo
                                 sod_shocktube_alecg_ale_sine.10.std.exo
                    BIN_RESULT out.e-s.0.1.0 out.e-s.10.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_sod_sine.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_diag.ndiff.cfg
                    LABELS alecg ale)

add_regression_test(compflow_euler_alecg_ale_sod_fluid ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES sod_shocktube_alecg_ale_fluid.q tube.exo
                    ARGS -c sod_shocktube_alecg_ale_fluid.q -v
                    BIN_BASELINE sod_shocktube_alecg_ale_fluid.0.std.exo
                                 sod_shocktube_alecg_ale_fluid.10.std.exo
                                 sod_shocktube_alecg_ale_fluid.20.std.exo
                    BIN_RESULT out.e-s.0.1.0 out.e-s.10.1.0 out.e-s.20.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_sod_fluid.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_diag.ndiff.cfg
                    LABELS alecg ale)

# Parallel

add_regression_test(compflow_euler_alecg_ale_gcl ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gcl.q rectangle_01_1.5k.exo
                    ARGS -c gcl.q -i rectangle_01_1.5k.exo -v
                          -u 0.5 +balancer RandCentLB +LBDebug 1
                    BIN_BASELINE gcl.0.std.exo.0
                                 gcl.0.std.exo.1
                                 gcl.0.std.exo.2
                                 gcl.0.std.exo.3
                                 gcl.0.std.exo.4
                                 gcl.0.std.exo.5
                                 gcl.0.std.exo.6
                                 gcl.10.std.exo.0
                                 gcl.10.std.exo.1
                                 gcl.10.std.exo.2
                                 gcl.10.std.exo.3
                                 gcl.10.std.exo.4
                                 gcl.10.std.exo.5
                                 gcl.10.std.exo.6
                    BIN_RESULT out.e-s.0.7.0
                               out.e-s.0.7.1
                               out.e-s.0.7.2
                               out.e-s.0.7.3
                               out.e-s.0.7.4
                               out.e-s.0.7.5
                               out.e-s.0.7.6
                               out.e-s.10.7.0
                               out.e-s.10.7.1
                               out.e-s.10.7.2
                               out.e-s.10.7.3
                               out.e-s.10.7.4
                               out.e-s.10.7.5
                               out.e-s.10.7.6
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_gcl.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gcl_diag.ndiff.cfg
                    LABELS alecg ale migration)

add_regression_test(compflow_euler_alecg_ale_sod_fluid ${INCITER_EXECUTABLE}
                    NUMPES 4 PPN 1
                    INPUTFILES sod_shocktube_alecg_ale_fluid.q tube.exo
                    ARGS -c sod_shocktube_alecg_ale_fluid.q -v
                          +balancer RandCentLB +LBDebug 1
                    BIN_BASELINE sod_shocktube_alecg_ale_fluid.0.std.exo.0
                                 sod_shocktube_alecg_ale_fluid.0.std.exo.1
                                 sod_shocktube_alecg_ale_fluid.0.std.exo.2
                                 sod_shocktube_alecg_ale_fluid.0.std.exo.3
                                 sod_shocktube_alecg_ale_fluid.10.std.exo.0
                                 sod_shocktube_alecg_ale_fluid.10.std.exo.1
                                 sod_shocktube_alecg_ale_fluid.10.std.exo.2
                                 sod_shocktube_alecg_ale_fluid.10.std.exo.3
                                 sod_shocktube_alecg_ale_fluid.20.std.exo.0
                                 sod_shocktube_alecg_ale_fluid.20.std.exo.1
                                 sod_shocktube_alecg_ale_fluid.20.std.exo.2
                                 sod_shocktube_alecg_ale_fluid.20.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                               out.e-s.10.4.0
                               out.e-s.10.4.1
                               out.e-s.10.4.2
                               out.e-s.10.4.3
                               out.e-s.20.4.0
                               out.e-s.20.4.1
                               out.e-s.20.4.2
                               out.e-s.20.4.3
                    BIN_DIFF_PROG_CONF exodiff_meshvel.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_sod_fluid.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF sod_diag.ndiff.cfg
                    LABELS alecg ale migration)
