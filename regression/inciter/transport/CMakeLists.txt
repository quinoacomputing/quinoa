# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(shear_centered_diffonly ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES shear_diffonly.q shear_centered_12k.exo
                               exodiff.cfg shear_centered_diffonly.std.exo
                    ARGS -c shear_diffonly.q -i shear_centered_12k.exo -v
                    BIN_BASELINE shear_centered_diffonly.std.exo
                    BIN_RESULT out.0
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(shear_centered_advdiff ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES shear_advdiff.q shear_centered_12k.exo
                               exodiff.cfg shear_centered_advdiff.std.exo
                    ARGS -c shear_advdiff.q -i shear_centered_12k.exo -v
                    BIN_BASELINE shear_centered_advdiff.std.exo
                    BIN_RESULT out.0
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(shear_centered_advdiffshear ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES shear_advdiffshear.q shear_centered_12k.exo
                               exodiff.cfg shear_centered_advdiffshear.std.exo
                    ARGS -c shear_advdiffshear.q -i shear_centered_12k.exo -v
                    BIN_BASELINE shear_centered_advdiffshear.std.exo
                    BIN_RESULT out.0
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(shear_centered_advdiffshear_c2 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES shear_advdiffshear_c2.q shear_centered_12k.exo
                               exodiff_c2.cfg
                               shear_centered_advdiffshear_c2.std.exo
                    ARGS -c shear_advdiffshear_c2.q -i shear_centered_12k.exo -v
                    BIN_BASELINE shear_centered_advdiffshear_c2.std.exo
                    BIN_RESULT out.0
                    BIN_DIFF_PROG_CONF exodiff_c2.cfg)


# Parallel + no virtualization

add_regression_test(shear_centered_advdiffshear ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES shear_advdiffshear.q shear_centered_12k.exo
                               exodiff.cfg
                               shear_centered_advdiffshear_pe4.std.exo.0
                               shear_centered_advdiffshear_pe4.std.exo.1
                               shear_centered_advdiffshear_pe4.std.exo.2
                               shear_centered_advdiffshear_pe4.std.exo.3
                    ARGS -c shear_advdiffshear.q -i shear_centered_12k.exo -v
                    BIN_BASELINE shear_centered_advdiffshear_pe4.std.exo.0
                                 shear_centered_advdiffshear_pe4.std.exo.1
                                 shear_centered_advdiffshear_pe4.std.exo.2
                                 shear_centered_advdiffshear_pe4.std.exo.3
                    BIN_RESULT out.0 out.1 out.2 out.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg)

# Parallel + virtualization

add_regression_test(shear_centered_advdiffshear_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES shear_advdiffshear.q shear_centered_12k.exo
                               exodiff_par.cfg
                               shear_centered_advdiffshear_pe4_u0.5.std.exo.0
                               shear_centered_advdiffshear_pe4_u0.5.std.exo.1
                               shear_centered_advdiffshear_pe4_u0.5.std.exo.2
                               shear_centered_advdiffshear_pe4_u0.5.std.exo.3
                               shear_centered_advdiffshear_pe4_u0.5.std.exo.4
                               shear_centered_advdiffshear_pe4_u0.5.std.exo.5
                               shear_centered_advdiffshear_pe4_u0.5.std.exo.6
                               shear_centered_advdiffshear_pe4_u0.5.std.exo.7
                    ARGS -c shear_advdiffshear.q -i shear_centered_12k.exo -v
                         -u 0.5
                    BIN_BASELINE shear_centered_advdiffshear_pe4_u0.5.std.exo.0
                                 shear_centered_advdiffshear_pe4_u0.5.std.exo.1
                                 shear_centered_advdiffshear_pe4_u0.5.std.exo.2
                                 shear_centered_advdiffshear_pe4_u0.5.std.exo.3
                                 shear_centered_advdiffshear_pe4_u0.5.std.exo.4
                                 shear_centered_advdiffshear_pe4_u0.5.std.exo.5
                                 shear_centered_advdiffshear_pe4_u0.5.std.exo.6
                                 shear_centered_advdiffshear_pe4_u0.5.std.exo.7
                    BIN_RESULT out.0 out.1 out.2 out.3 out.4 out.5 out.6 out.7
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff_par.cfg)
