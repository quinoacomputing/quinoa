# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(compflow_euler_taylor_green ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES taylor_green.q unitcube_1k.exo
                               exodiff.cfg taylor_green.std.exo
                    ARGS -c taylor_green.q -i unitcube_1k.exo -v
                    BIN_BASELINE taylor_green.std.exo
                    BIN_RESULT out.0
                    BIN_DIFF_PROG_CONF exodiff.cfg)

# Parallel + no virtualization

add_regression_test(compflow_euler_taylor_green ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES taylor_green.q unitcube_1k.exo
                               exodiff.cfg
                               taylor_green_pe4.std.exo.0
                               taylor_green_pe4.std.exo.1
                               taylor_green_pe4.std.exo.2
                               taylor_green_pe4.std.exo.3
                    ARGS -c taylor_green.q -i unitcube_1k.exo -v
                    BIN_BASELINE taylor_green_pe4.std.exo.0
                                 taylor_green_pe4.std.exo.1
                                 taylor_green_pe4.std.exo.2
                                 taylor_green_pe4.std.exo.3
                    BIN_RESULT out.0 out.1 out.2 out.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg)

# Parallel + virtualization

add_regression_test(compflow_euler_taylor_green_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES taylor_green.q unitcube_1k.exo
                               exodiff.cfg
                               taylor_green_pe4_u0.5.std.exo.0
                               taylor_green_pe4_u0.5.std.exo.1
                               taylor_green_pe4_u0.5.std.exo.2
                               taylor_green_pe4_u0.5.std.exo.3
                               taylor_green_pe4_u0.5.std.exo.4
                               taylor_green_pe4_u0.5.std.exo.5
                               taylor_green_pe4_u0.5.std.exo.6
                               taylor_green_pe4_u0.5.std.exo.7
                    ARGS -c taylor_green.q -i unitcube_1k.exo -v -u 0.5
                    BIN_BASELINE taylor_green_pe4_u0.5.std.exo.0
                                 taylor_green_pe4_u0.5.std.exo.1
                                 taylor_green_pe4_u0.5.std.exo.2
                                 taylor_green_pe4_u0.5.std.exo.3
                                 taylor_green_pe4_u0.5.std.exo.4
                                 taylor_green_pe4_u0.5.std.exo.5
                                 taylor_green_pe4_u0.5.std.exo.6
                                 taylor_green_pe4_u0.5.std.exo.7
                    BIN_RESULT out.0 out.1 out.2 out.3 out.4 out.5 out.6 out.7
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg)
