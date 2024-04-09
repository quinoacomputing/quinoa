# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# See regression/inciter/mesh_refinement/README.txt for naming conventions

# Serial

add_regression_test(amr_t0ref_u_trans_reord ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_amr_reord.lua unitsquare_01_955.exo
                    ARGS -c slot_cyl_amr_reord.lua
                         -i unitsquare_01_955.exo -v
                    BIN_BASELINE slot_cyl_pe1_u0.0.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF slot_cyl.ndiff.cfg
                    LABELS amr)

add_regression_test(amr_t0ref_u_compflow_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow_dg.lua unitcube_1k.exo
                    ARGS -c vortical_flow_dg.lua -i unitcube_1k.exo -v
                    BIN_BASELINE vortical_flow_dg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF vortical_flow_diag.ndiff.cfg
                    LABELS dg amr)

add_regression_test(amr_t0ref_u_compflow_reord_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES vortical_flow_reord_dg.lua unitcube_1k.exo
                    ARGS -c vortical_flow_reord_dg.lua -i unitcube_1k.exo -v
                    BIN_BASELINE vortical_flow_dg.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF vortical_flow_diag.ndiff.cfg
                    LABELS dg amr)

add_regression_test(amr_t0ref_cc_trans ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_amr_coords.lua unitcube_01_364.exo
                    ARGS -c slot_cyl_amr_coords.lua -i unitcube_01_364.exo -v
                    BIN_BASELINE amr_init_coords.1.std.exo.0
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS amr)

add_regression_test(amr_t0ref_u_trans_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES gauss_hump_dg.lua unitsquare_01_955_ss3.exo
                    ARGS -c gauss_hump_dg.lua -i unitsquare_01_955_ss3.exo -v
                    BIN_BASELINE gauss_hump_dg_t0ref.std.e-s.0.1.0
                                 gauss_hump_dg_t0ref.std.e-s.1.1.0
                                 gauss_hump_dg.std.exo
                    BIN_RESULT t0ref.e-s.0.1.0
                               t0ref.e-s.1.1.0
                               out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_gauss_hump_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
                    LABELS dg amr)

add_regression_test(amr_t0ref_u_trans_reord_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES gauss_hump_reord_dg.lua unitsquare_01_955_ss3.exo
                    ARGS -c gauss_hump_reord_dg.lua -i unitsquare_01_955_ss3.exo
                         -v
                    BIN_BASELINE gauss_hump_dg_t0ref.std.e-s.0.1.0
                                 gauss_hump_dg_t0ref.std.e-s.1.1.0
                                 gauss_hump_dg.std.exo
                    BIN_RESULT t0ref.e-s.0.1.0
                               t0ref.e-s.1.1.0
                               out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_gauss_hump_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
                    LABELS dg amr)

add_regression_test(amr_t0ref_udu_trans_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES gauss_hump_dg_uniform_deref.lua
                               unitsquare_01_955_ss3.exo
                    ARGS -c gauss_hump_dg_uniform_deref.lua
                         -i unitsquare_01_955_ss3.exo -v
                    BIN_BASELINE gauss_hump_dg_uniform_deref_t0ref.std.e-s.0.1.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.1.1.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.2.1.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.3.1.0
                                 gauss_hump_dg.std.exo
                    BIN_RESULT t0ref.e-s.0.1.0
                               t0ref.e-s.1.1.0
                               t0ref.e-s.2.1.0
                               t0ref.e-s.3.1.0
                               out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_gauss_hump_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
                    LABELS dg amr)

add_regression_test(amr_t0ref_ududu_trans_dg ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES gauss_hump_dg_uniform_deref_x2.lua
                               unitsquare_01_955_ss3.exo
                    ARGS -c gauss_hump_dg_uniform_deref_x2.lua
                         -i unitsquare_01_955_ss3.exo -v
                    BIN_BASELINE gauss_hump_dg_uniform_deref_t0ref.std.e-s.0.1.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.1.1.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.2.1.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.3.1.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.4.1.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.5.1.0
                                 gauss_hump_dg.std.exo
                    BIN_RESULT t0ref.e-s.0.1.0
                               t0ref.e-s.1.1.0
                               t0ref.e-s.2.1.0
                               t0ref.e-s.3.1.0
                               t0ref.e-s.4.1.0
                               t0ref.e-s.5.1.0
                               out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_gauss_hump_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
                    LABELS dg amr)

add_regression_test(amr_t0ref_idiu_trans ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES ic_uniform_deref.lua unitsquare_01_141.exo
                    ARGS -c ic_uniform_deref.lua -i unitsquare_01_141.exo -v
                    BIN_BASELINE ic_uniform_deref_t0ref.e-s.0.1.0
                                 ic_uniform_deref_t0ref.e-s.1.1.0
                                 ic_uniform_deref_t0ref.e-s.2.1.0
                                 ic_uniform_deref_t0ref.e-s.3.1.0
                                 ic_uniform_deref_t0ref.e-s.4.1.0
                                 ic_uniform_deref.std.exo
                    BIN_RESULT t0ref.e-s.0.1.0
                               t0ref.e-s.1.1.0
                               t0ref.e-s.2.1.0
                               t0ref.e-s.3.1.0
                               t0ref.e-s.4.1.0
                               out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    LABELS amr)

add_regression_test(amr_t0ref_uudd_trans ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_uudd.lua unitsquare_01_141.exo
                    ARGS -c slot_cyl_uudd.lua -i unitsquare_01_141.exo -v
                    BIN_BASELINE slot_cyl_uudd_t0ref.e-s.0.1.0
                                 slot_cyl_uudd_t0ref.e-s.1.1.0
                                 slot_cyl_uudd_t0ref.e-s.2.1.0
                                 slot_cyl_uudd_t0ref.e-s.3.1.0
                                 slot_cyl_uudd_t0ref.e-s.4.1.0
                                 slot_cyl.std.exo
                    BIN_RESULT t0ref.e-s.0.1.0
                               t0ref.e-s.1.1.0
                               t0ref.e-s.2.1.0
                               t0ref.e-s.3.1.0
                               t0ref.e-s.4.1.0
                               out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    LABELS amr)

# Parallel, no virtualization

add_regression_test(amr_t0ref_ii ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES ic_ic.lua unitsquare_01_141.exo
                    ARGS -c ic_ic.lua -i unitsquare_01_141.exo -v
                    BIN_BASELINE ic_ic_t0ref.e-s.0.2.0
                                 ic_ic_t0ref.e-s.0.2.1
                                 ic_ic_t0ref.e-s.1.2.0
                                 ic_ic_t0ref.e-s.1.2.1
                                 ic_ic_t0ref.e-s.2.2.0
                                 ic_ic_t0ref.e-s.2.2.1
                                 amr_ic_ic.std.exo.0
                                 amr_ic_ic.std.exo.1
                    BIN_RESULT t0ref.e-s.0.2.0
                               t0ref.e-s.0.2.1
                               t0ref.e-s.1.2.0
                               t0ref.e-s.1.2.1
                               t0ref.e-s.2.2.0
                               t0ref.e-s.2.2.1
                               out.e-s.0.2.0
                               out.e-s.0.2.1
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS amr)

add_regression_test(amr_t0ref_iu_trans ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES ic_uniform.lua unitsquare_01_141.exo
                    ARGS -c ic_uniform.lua -i unitsquare_01_141.exo -v
                    BIN_BASELINE amr_ic_uniform.std.exo.0
                                 amr_ic_uniform.std.exo.1
                    BIN_RESULT  out.e-s.0.2.0
                                out.e-s.0.2.1
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS amr)

add_regression_test(amr_t0ref_idiu_trans ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES ic_uniform_deref.lua unitsquare_01_141.exo
                    ARGS -c ic_uniform_deref.lua -i unitsquare_01_141.exo -v
                    BIN_BASELINE ic_uniform_deref_t0ref.e-s.0.2.0
                                 ic_uniform_deref_t0ref.e-s.0.2.1
                                 ic_uniform_deref_t0ref.e-s.1.2.0
                                 ic_uniform_deref_t0ref.e-s.1.2.1
                                 ic_uniform_deref_t0ref.e-s.2.2.0
                                 ic_uniform_deref_t0ref.e-s.2.2.1
                                 ic_uniform_deref_t0ref.e-s.3.2.0
                                 ic_uniform_deref_t0ref.e-s.3.2.1
                                 ic_uniform_deref_t0ref.e-s.4.2.0
                                 ic_uniform_deref_t0ref.e-s.4.2.1
                                 ic_uniform_deref.std.exo.0
                                 ic_uniform_deref.std.exo.1
                    BIN_RESULT t0ref.e-s.0.2.0
                               t0ref.e-s.0.2.1
                               t0ref.e-s.1.2.0
                               t0ref.e-s.1.2.1
                               t0ref.e-s.2.2.0
                               t0ref.e-s.2.2.1
                               t0ref.e-s.3.2.0
                               t0ref.e-s.3.2.1
                               t0ref.e-s.4.2.0
                               t0ref.e-s.4.2.1
                               out.e-s.0.2.0
                               out.e-s.0.2.1
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    LABELS amr)

add_regression_test(amr_t0ref_idddiu_trans ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES ic_uniform_deref2.lua unitsquare_01_141.exo
                    ARGS -c ic_uniform_deref2.lua -i unitsquare_01_141.exo -v
                    BIN_BASELINE ic_uniform_deref_t0ref.e-s.0.2.0
                                 ic_uniform_deref_t0ref.e-s.0.2.1
                                 ic_uniform_deref_t0ref.e-s.1.2.0
                                 ic_uniform_deref_t0ref.e-s.1.2.1
                                 ic_uniform_deref_t0ref.e-s.2.2.0
                                 ic_uniform_deref_t0ref.e-s.2.2.1
                                 ic_uniform_deref_t0ref.e-s.3.2.0
                                 ic_uniform_deref_t0ref.e-s.3.2.1
                                 ic_uniform_deref_t0ref.e-s.4.2.0
                                 ic_uniform_deref_t0ref.e-s.4.2.1
                                 ic_uniform_deref.std.exo.0
                                 ic_uniform_deref.std.exo.1
                    BIN_RESULT t0ref.e-s.0.2.0
                               t0ref.e-s.0.2.1
                               t0ref.e-s.1.2.0
                               t0ref.e-s.1.2.1
                               t0ref.e-s.4.2.0  # skip two no-op derefs
                               t0ref.e-s.4.2.1
                               t0ref.e-s.5.2.0
                               t0ref.e-s.5.2.1
                               t0ref.e-s.6.2.0
                               t0ref.e-s.6.2.1
                               out.e-s.0.2.0
                               out.e-s.0.2.1
                    BIN_DIFF_PROG_CONF exodiff_coords_but_ignore_time.cfg
                    BIN_DIFF_PROG_ARGS -m
                    LABELS amr)


add_regression_test(amr_t0ref_u_trans_dg ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gauss_hump_dg.lua unitsquare_01_955_ss3.exo
                    ARGS -c gauss_hump_dg.lua -i unitsquare_01_955_ss3.exo -v
                    BIN_BASELINE gauss_hump_dg_t0ref.std.e-s.0.4.0
                                 gauss_hump_dg_t0ref.std.e-s.0.4.1
                                 gauss_hump_dg_t0ref.std.e-s.0.4.2
                                 gauss_hump_dg_t0ref.std.e-s.0.4.3
                                 gauss_hump_dg_t0ref.std.e-s.1.4.0
                                 gauss_hump_dg_t0ref.std.e-s.1.4.1
                                 gauss_hump_dg_t0ref.std.e-s.1.4.2
                                 gauss_hump_dg_t0ref.std.e-s.1.4.3
                                 gauss_hump_dg_pe4.std.e-s.0.4.0
                                 gauss_hump_dg_pe4.std.e-s.0.4.1
                                 gauss_hump_dg_pe4.std.e-s.0.4.2
                                 gauss_hump_dg_pe4.std.e-s.0.4.3
                    BIN_RESULT t0ref.e-s.0.4.0
                               t0ref.e-s.0.4.1
                               t0ref.e-s.0.4.2
                               t0ref.e-s.0.4.3
                               t0ref.e-s.1.4.0
                               t0ref.e-s.1.4.1
                               t0ref.e-s.1.4.2
                               t0ref.e-s.1.4.3
                               out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.cfg
                                       exodiff_gauss_hump_dg.cfg
                                       exodiff_gauss_hump_dg.cfg
                                       exodiff_gauss_hump_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_gauss_hump_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
                    LABELS dg amr)

add_regression_test(amr_t0ref_udu_trans_dg ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gauss_hump_dg_uniform_deref.lua
                               unitsquare_01_955_ss3.exo
                    ARGS -c gauss_hump_dg_uniform_deref.lua
                         -i unitsquare_01_955_ss3.exo -v
                    BIN_BASELINE gauss_hump_dg_uniform_deref_t0ref.std.e-s.0.4.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.0.4.1
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.0.4.2
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.0.4.3
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.1.4.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.1.4.1
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.1.4.2
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.1.4.3
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.2.4.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.2.4.1
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.2.4.2
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.2.4.3
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.3.4.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.3.4.1
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.3.4.2
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.3.4.3
                                 gauss_hump_dg.std.exo.0
                                 gauss_hump_dg.std.exo.1
                                 gauss_hump_dg.std.exo.2
                                 gauss_hump_dg.std.exo.3
                    BIN_RESULT t0ref.e-s.0.4.0
                               t0ref.e-s.0.4.1
                               t0ref.e-s.0.4.2
                               t0ref.e-s.0.4.3
                               t0ref.e-s.1.4.0
                               t0ref.e-s.1.4.1
                               t0ref.e-s.1.4.2
                               t0ref.e-s.1.4.3
                               t0ref.e-s.2.4.0
                               t0ref.e-s.2.4.1
                               t0ref.e-s.2.4.2
                               t0ref.e-s.2.4.3
                               t0ref.e-s.3.4.0
                               t0ref.e-s.3.4.1
                               t0ref.e-s.3.4.2
                               t0ref.e-s.3.4.3
                               out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.cfg
                                       exodiff_gauss_hump_dg.cfg
                                       exodiff_gauss_hump_dg.cfg
                                       exodiff_gauss_hump_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_gauss_hump_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
                    LABELS dg amr)

add_regression_test(amr_t0ref_ududu_trans_dg ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES gauss_hump_dg_uniform_deref_x2.lua
                               unitsquare_01_955_ss3.exo
                    ARGS -c gauss_hump_dg_uniform_deref_x2.lua
                         -i unitsquare_01_955_ss3.exo -v
                    BIN_BASELINE gauss_hump_dg_uniform_deref_t0ref.std.e-s.0.4.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.0.4.1
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.0.4.2
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.0.4.3
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.1.4.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.1.4.1
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.1.4.2
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.1.4.3
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.2.4.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.2.4.1
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.2.4.2
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.2.4.3
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.3.4.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.3.4.1
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.3.4.2
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.3.4.3
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.4.4.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.4.4.1
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.4.4.2
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.4.4.3
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.5.4.0
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.5.4.1
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.5.4.2
                                 gauss_hump_dg_uniform_deref_t0ref.std.e-s.5.4.3
                                 gauss_hump_dg.std.exo.0
                                 gauss_hump_dg.std.exo.1
                                 gauss_hump_dg.std.exo.2
                                 gauss_hump_dg.std.exo.3
                    BIN_RESULT t0ref.e-s.0.4.0
                               t0ref.e-s.0.4.1
                               t0ref.e-s.0.4.2
                               t0ref.e-s.0.4.3
                               t0ref.e-s.1.4.0
                               t0ref.e-s.1.4.1
                               t0ref.e-s.1.4.2
                               t0ref.e-s.1.4.3
                               t0ref.e-s.2.4.0
                               t0ref.e-s.2.4.1
                               t0ref.e-s.2.4.2
                               t0ref.e-s.2.4.3
                               t0ref.e-s.3.4.0
                               t0ref.e-s.3.4.1
                               t0ref.e-s.3.4.2
                               t0ref.e-s.3.4.3
                               t0ref.e-s.4.4.0
                               t0ref.e-s.4.4.1
                               t0ref.e-s.4.4.2
                               t0ref.e-s.4.4.3
                               t0ref.e-s.5.4.0
                               t0ref.e-s.5.4.1
                               t0ref.e-s.5.4.2
                               t0ref.e-s.5.4.3
                               out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_CONF exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.t0ref.cfg
                                       exodiff_gauss_hump_dg.cfg
                                       exodiff_gauss_hump_dg.cfg
                                       exodiff_gauss_hump_dg.cfg
                                       exodiff_gauss_hump_dg.cfg
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_gauss_hump_dg.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF gauss_hump_diag.ndiff.cfg
                    LABELS dg amr)

add_regression_test(amr_t0ref_uudd_trans ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES slot_cyl_uudd.lua unitsquare_01_141.exo
                    ARGS -c slot_cyl_uudd.lua -i unitsquare_01_141.exo -v
                    BIN_BASELINE slot_cyl_uudd_t0ref.e-s.0.2.0
                                 slot_cyl_uudd_t0ref.e-s.0.2.1
                                 slot_cyl_uudd_t0ref.e-s.1.2.0
                                 slot_cyl_uudd_t0ref.e-s.1.2.1
                                 slot_cyl_uudd_t0ref.e-s.2.2.0
                                 slot_cyl_uudd_t0ref.e-s.2.2.1
                                 slot_cyl_uudd_t0ref.e-s.3.2.0
                                 slot_cyl_uudd_t0ref.e-s.3.2.1
                                 slot_cyl_uudd_t0ref.e-s.4.2.0
                                 slot_cyl_uudd_t0ref.e-s.4.2.1
                                 slot_cyl.std.exo.0
                                 slot_cyl.std.exo.1
                    BIN_RESULT t0ref.e-s.0.2.0
                               t0ref.e-s.0.2.1
                               t0ref.e-s.1.2.0
                               t0ref.e-s.1.2.1
                               t0ref.e-s.2.2.0
                               t0ref.e-s.2.2.1
                               t0ref.e-s.3.2.0
                               t0ref.e-s.3.2.1
                               t0ref.e-s.4.2.0
                               t0ref.e-s.4.2.1
                               out.e-s.0.2.0
                               out.e-s.0.2.1
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    LABELS amr)

add_regression_test(amr_t0ref_uuuddd_trans ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES slot_cyl_uuuddd.lua unitsquare_01_141.exo
                    ARGS -c slot_cyl_uuuddd.lua -i unitsquare_01_141.exo -v
                    BIN_BASELINE slot_cyl_uuuddd_t0ref.e-s.0.2.0
                                 slot_cyl_uuuddd_t0ref.e-s.0.2.1
                                 slot_cyl_uuuddd_t0ref.e-s.1.2.0
                                 slot_cyl_uuuddd_t0ref.e-s.1.2.1
                                 slot_cyl_uuuddd_t0ref.e-s.2.2.0
                                 slot_cyl_uuuddd_t0ref.e-s.2.2.1
                                 slot_cyl_uuuddd_t0ref.e-s.3.2.0
                                 slot_cyl_uuuddd_t0ref.e-s.3.2.1
                                 slot_cyl_uuuddd_t0ref.e-s.4.2.0
                                 slot_cyl_uuuddd_t0ref.e-s.4.2.1
                                 slot_cyl_uuuddd_t0ref.e-s.5.2.0
                                 slot_cyl_uuuddd_t0ref.e-s.5.2.1
                                 slot_cyl_uuuddd_t0ref.e-s.6.2.0
                                 slot_cyl_uuuddd_t0ref.e-s.6.2.1
                                 slot_cyl.std.exo.0
                                 slot_cyl.std.exo.1
                    BIN_RESULT t0ref.e-s.0.2.0
                               t0ref.e-s.0.2.1
                               t0ref.e-s.1.2.0
                               t0ref.e-s.1.2.1
                               t0ref.e-s.2.2.0
                               t0ref.e-s.2.2.1
                               t0ref.e-s.3.2.0
                               t0ref.e-s.3.2.1
                               t0ref.e-s.4.2.0
                               t0ref.e-s.4.2.1
                               t0ref.e-s.5.2.0
                               t0ref.e-s.5.2.1
                               t0ref.e-s.6.2.0
                               t0ref.e-s.6.2.1
                               out.e-s.0.2.0
                               out.e-s.0.2.1
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    BIN_DIFF_PROG_ARGS -m
                    LABELS amr)
