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
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gcl_diag.ndiff.cfg
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
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gcl_diag.ndiff.cfg
                    LABELS alecg ale migration)
