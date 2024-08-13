# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(slotted_cylinder_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_dg.lua unitcube_01_31k.exo
                    ARGS -c slot_cyl_dg.lua -i unitcube_01_31k.exo -v
                    BIN_BASELINE slot_cyl_dg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    TEXT_BASELINE diag_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF slo_cyl_diag.ndiff.cfg
                    LABELS dg)

#add_regression_test(slotted_cylinder_dgp1 ${INCITER_EXECUTABLE}
#                    NUMPES 1
#                    INPUTFILES slot_cyl_dgp1.lua unitcube_01_31k.exo
#                    ARGS -c slot_cyl_dgp1.lua -i unitcube_01_31k.exo -v
#                    BIN_BASELINE slot_cyl_dgp1.std.exo
#                    BIN_RESULT out.e-s.0.1.0
#                    BIN_DIFF_PROG_ARGS -m
#                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
#                    TEXT_BASELINE diag_dgp1.std
#                    TEXT_RESULT diag
#                    TEXT_DIFF_PROG_CONF slo_cly_diag.ndiff.cfg
#                    LABELS dg)

# Parallel

add_regression_test(slotted_cylinder_dg ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl_dg.lua unitcube_01_31k.exo
                    ARGS -c slot_cyl_dg.lua -i unitcube_01_31k.exo -v
                    BIN_BASELINE slot_cyl_dg_pe4.std.exo.0
                                 slot_cyl_dg_pe4.std.exo.1
                                 slot_cyl_dg_pe4.std.exo.2
                                 slot_cyl_dg_pe4.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    TEXT_BASELINE diag_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF slo_cyl_diag.ndiff.cfg
                    LABELS dg)

#add_regression_test(slotted_cylinder_dgp1 ${INCITER_EXECUTABLE}
#                    NUMPES 4
#                    INPUTFILES slot_cyl_dgp1.lua unitcube_01_31k.exo
#                    ARGS -c slot_cyl_dgp1.lua -i unitcube_01_31k.exo -v
#                    BIN_BASELINE slot_cyl_dgp1_pe4.std.exo.0
#                                 slot_cyl_dgp1_pe4.std.exo.1
#                                 slot_cyl_dgp1_pe4.std.exo.2
#                                 slot_cyl_dgp1_pe4.std.exo.3
#                    BIN_RESULT out.e-s.0.4.0
#                               out.e-s.0.4.1
#                               out.e-s.0.4.2
#                               out.e-s.0.4.3
#                    BIN_DIFF_PROG_ARGS -m
#                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
#                    TEXT_BASELINE diag_dgp1.std
#                    TEXT_RESULT diag
#                    TEXT_DIFF_PROG_CONF slo_cyl_diag.ndiff.cfg
#                    LABELS dg)
