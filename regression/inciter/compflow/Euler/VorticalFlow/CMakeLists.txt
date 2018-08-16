# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial

add_regression_test(compflow_euler_vorticalflow ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow.q unitcube_1k.exo diag.std
                               exodiff.cfg vortical_flow.std.exo
                    ARGS -c vortical_flow.q -i unitcube_1k.exo -v
                    BIN_BASELINE vortical_flow.std.exo
                    BIN_RESULT out.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF vortical_flow_diag.ndiff.cfg)

add_regression_test(compflow_euler_vorticalflow_diagcg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow_diagcg.q unitcube_1k.exo
                               diag_diagcg.std
                               exodiff.cfg vortical_flow_diagcg.std.exo
                    ARGS -c vortical_flow_diagcg.q -i unitcube_1k.exo -v
                    BIN_BASELINE vortical_flow_diagcg.std.exo
                    BIN_RESULT out.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_diagcg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF vortical_flow_diag.ndiff.cfg)

add_regression_test(compflow_euler_vorticalflow_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow_dg.q unitcube_1k.exo
                               diag_dg.std
                               exodiff_dg.cfg vortical_flow_dg.std.exo
                    ARGS -c vortical_flow_dg.q -i unitcube_1k.exo -v
                    BIN_BASELINE vortical_flow_dg.std.exo
                    BIN_RESULT out.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF vortical_flow_diag.ndiff.cfg)

add_regression_test(compflow_euler_vorticalflow_dg_lf ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow_dg_lf.q unitcube_1k.exo
                               diag_dg_lf.std
                               exodiff_dg.cfg vortical_flow_dg_lf.std.exo
                    ARGS -c vortical_flow_dg_lf.q -i unitcube_1k.exo -v
                    BIN_BASELINE vortical_flow_dg_lf.std.exo
                    BIN_RESULT out.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dg_lf.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF vortical_flow_diag.ndiff.cfg)

# Parallel + no virtualization

add_regression_test(compflow_euler_vorticalflow ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES vortical_flow.q unitcube_1k.exo
                               exodiff.cfg
                               diag.std
                               vortical_flow_pe4.std.exo.0
                               vortical_flow_pe4.std.exo.1
                               vortical_flow_pe4.std.exo.2
                               vortical_flow_pe4.std.exo.3
                    ARGS -c vortical_flow.q -i unitcube_1k.exo -v
                    BIN_BASELINE vortical_flow_pe4.std.exo.0
                                 vortical_flow_pe4.std.exo.1
                                 vortical_flow_pe4.std.exo.2
                                 vortical_flow_pe4.std.exo.3
                    BIN_RESULT out.0 out.1 out.2 out.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF vortical_flow_diag.ndiff.cfg)

# Parallel + virtualization

add_regression_test(compflow_euler_vorticalflow_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES vortical_flow.q unitcube_1k.exo
                               exodiff.cfg
                               diag.std
                               vortical_flow_pe4_u0.5.std.exo.0
                               vortical_flow_pe4_u0.5.std.exo.1
                               vortical_flow_pe4_u0.5.std.exo.2
                               vortical_flow_pe4_u0.5.std.exo.3
                               vortical_flow_pe4_u0.5.std.exo.4
                               vortical_flow_pe4_u0.5.std.exo.5
                               vortical_flow_pe4_u0.5.std.exo.6
                               vortical_flow_pe4_u0.5.std.exo.7
                    ARGS -c vortical_flow.q -i unitcube_1k.exo -v -u 0.5
                    BIN_BASELINE vortical_flow_pe4_u0.5.std.exo.0
                                 vortical_flow_pe4_u0.5.std.exo.1
                                 vortical_flow_pe4_u0.5.std.exo.2
                                 vortical_flow_pe4_u0.5.std.exo.3
                                 vortical_flow_pe4_u0.5.std.exo.4
                                 vortical_flow_pe4_u0.5.std.exo.5
                                 vortical_flow_pe4_u0.5.std.exo.6
                                 vortical_flow_pe4_u0.5.std.exo.7
                    BIN_RESULT out.0 out.1 out.2 out.3 out.4 out.5 out.6 out.7
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF vortical_flow_diag.ndiff.cfg)

add_regression_test(compflow_euler_vorticalflow_diagcg_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES vortical_flow_diagcg.q unitcube_1k.exo
                               diag_diagcg.std exodiff.cfg
                               vortical_flow_diagcg_pe4_u0.5.std.exo.0
                               vortical_flow_diagcg_pe4_u0.5.std.exo.1
                               vortical_flow_diagcg_pe4_u0.5.std.exo.2
                               vortical_flow_diagcg_pe4_u0.5.std.exo.3
                               vortical_flow_diagcg_pe4_u0.5.std.exo.4
                               vortical_flow_diagcg_pe4_u0.5.std.exo.5
                               vortical_flow_diagcg_pe4_u0.5.std.exo.6
                               vortical_flow_diagcg_pe4_u0.5.std.exo.7
                    ARGS -c vortical_flow_diagcg.q -i unitcube_1k.exo -v -u 0.5
                    BIN_BASELINE vortical_flow_diagcg_pe4_u0.5.std.exo.0
                                 vortical_flow_diagcg_pe4_u0.5.std.exo.1
                                 vortical_flow_diagcg_pe4_u0.5.std.exo.2
                                 vortical_flow_diagcg_pe4_u0.5.std.exo.3
                                 vortical_flow_diagcg_pe4_u0.5.std.exo.4
                                 vortical_flow_diagcg_pe4_u0.5.std.exo.5
                                 vortical_flow_diagcg_pe4_u0.5.std.exo.6
                                 vortical_flow_diagcg_pe4_u0.5.std.exo.7
                    BIN_RESULT out.0 out.1 out.2 out.3 out.4 out.5 out.6 out.7
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_diagcg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF vortical_flow_diag.ndiff.cfg)
