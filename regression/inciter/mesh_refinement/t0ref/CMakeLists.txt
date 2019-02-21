# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# See regression/inciter/mesh_refinement/README.txt for naming conventions

# Virtualization = 0

add_regression_test(amr_t0ref_u_trans_diagcg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_amr_diagcg.q unitsquare_01_955.exo
                    ARGS -c slot_cyl_amr_diagcg.q -i unitsquare_01_955.exo -v
                    BIN_BASELINE slot_cyl_diagcg_pe1_u0.0.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_diagcg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF slot_cyl_diagcg.ndiff.cfg
                    LABELS amr)

add_regression_test(amr_t0ref_u_trans_reord_diagcg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_amr_reord_diagcg.q unitsquare_01_955.exo
                    ARGS -c slot_cyl_amr_reord_diagcg.q
                         -i unitsquare_01_955.exo -v
                    BIN_BASELINE slot_cyl_diagcg_pe1_u0.0.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_diagcg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF slot_cyl_diagcg.ndiff.cfg
                    LABELS amr)

add_regression_test(amr_t0ref_u_compflow_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow_dg.q unitcube_1k.exo
                    ARGS -c vortical_flow_dg.q -i unitcube_1k.exo -v
                    BIN_BASELINE vortical_flow_dg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF vortical_flow_diag.ndiff.cfg
                    LABELS amr)

add_regression_test(amr_t0ref_u_compflow_reord_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow_reord_dg.q unitcube_1k.exo
                    ARGS -c vortical_flow_reord_dg.q -i unitcube_1k.exo -v
                    BIN_BASELINE vortical_flow_dg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF vortical_flow_diag.ndiff.cfg
                    LABELS amr)

add_regression_test(amr_t0ref_cc ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_amr_coords.q unitcube_01_364.exo
                    ARGS -c slot_cyl_amr_coords.q -i unitcube_01_364.exo -v
                    BIN_BASELINE amr_init_coords.1.std.exo.0
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS amr)

# Parallel, no virtualization

add_regression_test(amr_t0ref_u_trans_diagcg ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl_amr_diagcg.q unitsquare_01_955.exo
                    ARGS -c slot_cyl_amr_diagcg.q -i unitsquare_01_955.exo -v
                    BIN_BASELINE slot_cyl_diagcg_pe4_u0.0.std.exo.0
                                 slot_cyl_diagcg_pe4_u0.0.std.exo.1
                                 slot_cyl_diagcg_pe4_u0.0.std.exo.2
                                 slot_cyl_diagcg_pe4_u0.0.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS amr)

add_regression_test(amr_t0ref_ii ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES ic_ic.q unitsquare_01_141.exo
                    ARGS -c ic_ic.q -i unitsquare_01_141.exo -v
                    BIN_BASELINE amr_ic_ic.std.exo.0
                                 amr_ic_ic.std.exo.1
                    # TODO: double check this numbering scheme
                    BIN_RESULT  out.e-s.0.2.0
                                out.e-s.0.2.1
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS amr)

add_regression_test(amr_t0ref_iu ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES ic_uniform.q unitsquare_01_141.exo
                    ARGS -c ic_uniform.q -i unitsquare_01_141.exo -v
                    BIN_BASELINE amr_ic_uniform.std.exo.0
                                 amr_ic_uniform.std.exo.1
                    # TODO: double check this numbering scheme
                    BIN_RESULT  out.e-s.0.2.0
                                out.e-s.0.2.1
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS amr)

# Parallel, virtualization = 0.5

add_regression_test(amr_t0ref_u_trans_diagcg_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl_amr_diagcg.q unitsquare_01_955.exo
                    ARGS -c slot_cyl_amr_diagcg.q -i unitsquare_01_955.exo
                         -v -u 0.5
                    BIN_BASELINE slot_cyl_diagcg_pe4_u0.5.std.exo.0
                                 slot_cyl_diagcg_pe4_u0.5.std.exo.1
                                 slot_cyl_diagcg_pe4_u0.5.std.exo.2
                                 slot_cyl_diagcg_pe4_u0.5.std.exo.3
                                 slot_cyl_diagcg_pe4_u0.5.std.exo.4
                                 slot_cyl_diagcg_pe4_u0.5.std.exo.5
                                 slot_cyl_diagcg_pe4_u0.5.std.exo.6
                                 slot_cyl_diagcg_pe4_u0.5.std.exo.7
                    BIN_RESULT out.e-s.0.8.0
                               out.e-s.0.8.1
                               out.e-s.0.8.2
                               out.e-s.0.8.3
                               out.e-s.0.8.4
                               out.e-s.0.8.5
                               out.e-s.0.8.6
                               out.e-s.0.8.7
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS amr)
