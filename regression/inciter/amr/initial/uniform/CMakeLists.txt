# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Virtualization = 0

add_regression_test(initial_uniform_amr ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_amr.q unitsquare_01_955.exo
                               exodiff.cfg slot_cyl_pe1_u0.0.std.exo
                    ARGS -c slot_cyl_amr.q -i unitsquare_01_955.exo -v
                    BIN_BASELINE slot_cyl_pe1_u0.0.std.exo
                    BIN_RESULT out.0
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg)

# Parallel, no virtualization

add_regression_test(initial_uniform_amr ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl_amr.q unitsquare_01_955.exo
                               exodiff.cfg
                               slot_cyl_pe4_u0.0.std.exo.0
                               slot_cyl_pe4_u0.0.std.exo.1
                               slot_cyl_pe4_u0.0.std.exo.2
                               slot_cyl_pe4_u0.0.std.exo.3
                    ARGS -c slot_cyl_amr.q -i unitsquare_01_955.exo -v
                    BIN_BASELINE slot_cyl_pe4_u0.0.std.exo.0
                                 slot_cyl_pe4_u0.0.std.exo.1
                                 slot_cyl_pe4_u0.0.std.exo.2
                                 slot_cyl_pe4_u0.0.std.exo.3
                    BIN_RESULT out.0
                               out.1
                               out.2
                               out.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg)

# Parallel, virtualization = 0.7

add_regression_test(initial_uniform_amr_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl_amr.q unitsquare_01_955.exo
                               exodiff.cfg
                               slot_cyl_pe4_u0.5.std.exo.0
                               slot_cyl_pe4_u0.5.std.exo.1
                               slot_cyl_pe4_u0.5.std.exo.2
                               slot_cyl_pe4_u0.5.std.exo.3
                               slot_cyl_pe4_u0.5.std.exo.4
                               slot_cyl_pe4_u0.5.std.exo.5
                               slot_cyl_pe4_u0.5.std.exo.6
                               slot_cyl_pe4_u0.5.std.exo.7
                    ARGS -c slot_cyl_amr.q -i unitsquare_01_955.exo -v -u 0.5
                    BIN_BASELINE slot_cyl_pe4_u0.5.std.exo.0
                                 slot_cyl_pe4_u0.5.std.exo.1
                                 slot_cyl_pe4_u0.5.std.exo.2
                                 slot_cyl_pe4_u0.5.std.exo.3
                                 slot_cyl_pe4_u0.5.std.exo.4
                                 slot_cyl_pe4_u0.5.std.exo.5
                                 slot_cyl_pe4_u0.5.std.exo.6
                                 slot_cyl_pe4_u0.5.std.exo.7
                    BIN_RESULT out.0
                               out.1
                               out.2
                               out.3
                               out.4
                               out.5
                               out.6
                               out.7
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg)
