# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(compflow_euler_nleg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES nleg.q unitcube_1k.exo
                               exodiff.cfg nleg.std.exo
                    ARGS -c nleg.q -i unitcube_1k.exo -v
                    BIN_BASELINE nleg.std.exo
                    BIN_RESULT out.0
                    BIN_DIFF_PROG_CONF exodiff.cfg)

# Parallel + no virtualization

add_regression_test(compflow_euler_nleg ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES nleg.q unitcube_1k.exo
                               exodiff.cfg
                               nleg_pe4.std.exo.0
                               nleg_pe4.std.exo.1
                               nleg_pe4.std.exo.2
                               nleg_pe4.std.exo.3
                    ARGS -c nleg.q -i unitcube_1k.exo -v
                    BIN_BASELINE nleg_pe4.std.exo.0
                                 nleg_pe4.std.exo.1
                                 nleg_pe4.std.exo.2
                                 nleg_pe4.std.exo.3
                    BIN_RESULT out.0 out.1 out.2 out.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg)

# Parallel + virtualization

add_regression_test(compflow_euler_nleg_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES nleg.q unitcube_1k.exo
                               exodiff.cfg
                               nleg_pe4_u0.5.std.exo.0
                               nleg_pe4_u0.5.std.exo.1
                               nleg_pe4_u0.5.std.exo.2
                               nleg_pe4_u0.5.std.exo.3
                               nleg_pe4_u0.5.std.exo.4
                               nleg_pe4_u0.5.std.exo.5
                               nleg_pe4_u0.5.std.exo.6
                               nleg_pe4_u0.5.std.exo.7
                    ARGS -c nleg.q -i unitcube_1k.exo -v -u 0.5
                    BIN_BASELINE nleg_pe4_u0.5.std.exo.0
                                 nleg_pe4_u0.5.std.exo.1
                                 nleg_pe4_u0.5.std.exo.2
                                 nleg_pe4_u0.5.std.exo.3
                                 nleg_pe4_u0.5.std.exo.4
                                 nleg_pe4_u0.5.std.exo.5
                                 nleg_pe4_u0.5.std.exo.6
                                 nleg_pe4_u0.5.std.exo.7
                    BIN_RESULT out.0 out.1 out.2 out.3 out.4 out.5 out.6 out.7
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg)
