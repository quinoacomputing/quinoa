# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Virtualization = 0

add_regression_test(fct ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v
                    BIN_BASELINE slot_cyl_pe1_u0.0.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS diagcg)

add_regression_test(fct ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -f
                    BIN_BASELINE slot_cyl_pe4_u0.0.std.exo.0
                                 slot_cyl_pe4_u0.0.std.exo.1
                                 slot_cyl_pe4_u0.0.std.exo.2
                                 slot_cyl_pe4_u0.0.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS diagcg)

# Virtualization = 0.5

add_regression_test(fct_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.5
                    BIN_BASELINE slot_cyl_pe1_u0.5.std.exo.0
                                 slot_cyl_pe1_u0.5.std.exo.1
                    BIN_RESULT out.e-s.0.2.0
                               out.e-s.0.2.1
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS diagcg)

# Virtualization = 0.9

add_regression_test(fct_u0.9 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.9
                    BIN_BASELINE slot_cyl_pe1_u0.9.std.exo.0
                                 slot_cyl_pe1_u0.9.std.exo.1
                                 slot_cyl_pe1_u0.9.std.exo.2
                                 slot_cyl_pe1_u0.9.std.exo.3
                                 slot_cyl_pe1_u0.9.std.exo.4
                                 slot_cyl_pe1_u0.9.std.exo.5
                                 slot_cyl_pe1_u0.9.std.exo.6
                                 slot_cyl_pe1_u0.9.std.exo.7
                                 slot_cyl_pe1_u0.9.std.exo.8
                    BIN_RESULT out.e-s.0.9.0
                               out.e-s.0.9.1
                               out.e-s.0.9.2
                               out.e-s.0.9.3
                               out.e-s.0.9.4
                               out.e-s.0.9.5
                               out.e-s.0.9.6
                               out.e-s.0.9.7
                               out.e-s.0.9.8
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS diagcg)
