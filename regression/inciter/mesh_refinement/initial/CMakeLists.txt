# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Naming convention for AMR regression tests
#
# Example: amr_init_u1_trans_reord_dg_u0.4
#
# Legend:
# * amr - exercising adaptive mesh refinement (AMR)
# * init - AMR before t=0 (initial AMR),
#          other possibilities instead of 'init':
#          - dt - AMR durging time stepping
# * u1 - uniform refinement (1:8 for all tets), single step
#        other possibilities instead of 'u':
#        - i - initial conditions based adaptive refinement,
#        - c - coordinate based refinement,
#        - l - refine a list of tagged edges
#        - C - combinations of the above
# * trans - type of physics, e.g., compflow, transport
#           other possibilities instead of 'trans':
#           - compflow - compressible flow
# * reord - perform PE-locality mesh node reordering during setup
# * dg - discontinuous Galerkin discretization scheme, 1st order, DG = DG(P0),
#        other possibilities instead of 'dg':
#        - dgp1 - DG(P1), 2nd-order DG
#        - diagcg - lumped-mass (diagonal) continuous Galerkin
# * u0.4 - Non-zero virtualization u = 0.4

# Virtualization = 0

add_regression_test(amr_init_u1_trans_diagcg ${INCITER_EXECUTABLE}
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

add_regression_test(amr_init_u1_trans_reord_diagcg ${INCITER_EXECUTABLE}
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

add_regression_test(amr_init_u1_compflow_dg ${INCITER_EXECUTABLE}
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

add_regression_test(amr_init_u1_compflow_reord_dg ${INCITER_EXECUTABLE}
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

# Parallel, no virtualization

add_regression_test(amr_init_u1_trans_diagcg ${INCITER_EXECUTABLE}
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

# Parallel, virtualization = 0.5

add_regression_test(amr_init_u1_trans_diagcg_u0.5 ${INCITER_EXECUTABLE}
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

add_regression_test(amr_init_coords ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_amr_coords.q unitcube_01_364.exo
                    ARGS -c slot_cyl_amr_coords.q -i unitcube_01_364.exo -v
                    BIN_BASELINE amr_init_coords.1.std.exo.0
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS amr)

add_regression_test(amr_init_ic_ic ${INCITER_EXECUTABLE}
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

add_regression_test(amr_init_ic_uniform ${INCITER_EXECUTABLE}
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
