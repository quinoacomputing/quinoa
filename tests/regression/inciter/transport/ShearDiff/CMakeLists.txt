# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(shear_centered_diffonly ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES shear_diffonly.q shear_centered_12k.exo
                    ARGS -c shear_diffonly.q -i shear_centered_12k.exo -v
                    BIN_BASELINE shear_centered_diffonly.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE shear_centered_diffonly.diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF shear.ndiff.cfg
                    LABELS diagcg)

add_regression_test(shear_centered_diffonly_nofct ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES shear_diffonly_nofct.q shear_centered_12k.exo
                    ARGS -c shear_diffonly_nofct.q -i shear_centered_12k.exo -v
                    BIN_BASELINE shear_centered_diffonly_nofct.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE shear_centered_diffonly_nofct.diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF shear.ndiff.cfg
                    LABELS diagcg)

add_regression_test(shear_centered_advdiff ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES shear_advdiff.q shear_centered_12k.exo
                    ARGS -c shear_advdiff.q -i shear_centered_12k.exo -v
                    BIN_BASELINE shear_centered_advdiff.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS diagcg)

add_regression_test(shear_centered_advdiffshear ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES shear_advdiffshear.q shear_centered_12k.exo
                    ARGS -c shear_advdiffshear.q -i shear_centered_12k.exo -v
                    BIN_BASELINE shear_centered_advdiffshear.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE shear_centered_advdiffshear.diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF shear.ndiff.cfg
                    LABELS diagcg)

add_regression_test(shear_centered_advdiffshear_c2 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES shear_advdiffshear_c2.q shear_centered_12k.exo
                    ARGS -c shear_advdiffshear_c2.q -i shear_centered_12k.exo -v
                    BIN_BASELINE shear_centered_advdiffshear_c2.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_c2.cfg
                    TEXT_BASELINE shear_centered_advdiffshear_c2.diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF shear.ndiff.cfg
                    LABELS diagcg)

# Parallel + no virtualization

add_regression_test(shear_centered_advdiffshear ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES shear_advdiffshear.q shear_centered_12k.exo
                    ARGS -c shear_advdiffshear.q -i shear_centered_12k.exo -v
                    BIN_BASELINE shear_centered_advdiffshear_pe4.std.exo.0
                                 shear_centered_advdiffshear_pe4.std.exo.1
                                 shear_centered_advdiffshear_pe4.std.exo.2
                                 shear_centered_advdiffshear_pe4.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff_par.cfg
                    TEXT_BASELINE shear_centered_advdiffshear.diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF shear.ndiff.cfg
                    LABELS diagcg)

# Parallel + virtualization

add_regression_test(shear_centered_advdiffshear_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES shear_advdiffshear.q shear_centered_12k.exo
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
                    BIN_RESULT out.e-s.0.8.0
                               out.e-s.0.8.1
                               out.e-s.0.8.2
                               out.e-s.0.8.3
                               out.e-s.0.8.4
                               out.e-s.0.8.5
                               out.e-s.0.8.6
                               out.e-s.0.8.7
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff_par.cfg
                    TEXT_BASELINE shear_centered_advdiffshear.diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF shear.ndiff.cfg
                    LABELS diagcg)
