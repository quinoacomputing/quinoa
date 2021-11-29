# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Virtualization = 0

add_regression_test(cfl ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_cfl.q unitsquare_01_3.6k.exo
                    ARGS -c slot_cyl_cfl.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE slot_cyl_cfl_pe1_u0.0.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS diagcg)

if (HAS_OMEGA_H)
  add_regression_test(cfl_osh ${INCITER_EXECUTABLE}
                      NUMPES 1
                      INPUTFILES slot_cyl_cfl.q box10.osh
                      ARGS -c slot_cyl_cfl.q -i box10.osh -v
                      BIN_BASELINE slot_cyl_cfl_osh10_pe1.std.exo
                      BIN_RESULT out.e-s.0.1.0
                      BIN_DIFF_PROG_CONF exodiff.cfg
                      LABELS diagcg)

  # parallel with nranks = nparts
  add_regression_test(cfl_osh ${INCITER_EXECUTABLE}
                      NUMPES 2 PPN 1
                      INPUTFILES slot_cyl_cfl.q box10_2.osh
                      ARGS -c slot_cyl_cfl.q -i box10_2.osh -v
                      BIN_BASELINE slot_cyl_cfl_osh10_pe2.std.exo.0
                                   slot_cyl_cfl_osh10_pe2.std.exo.1
                      BIN_RESULT out.e-s.0.2.0 out.e-s.0.2.1
                      BIN_DIFF_PROG_CONF exodiff.cfg
                      LABELS diagcg)

#  # parallel with nranks > nparts
#  add_regression_test(cfl_osh ${INCITER_EXECUTABLE}
#                      NUMPES 4
#                      INPUTFILES slot_cyl_cfl.q box10_2.osh
#                      ARGS -c slot_cyl_cfl.q -i box10_2.osh -v
#                      BIN_BASELINE slot_cyl_cfl_osh10_pe4.std.exo.0
#                                   slot_cyl_cfl_osh10_pe4.std.exo.1
#                                   slot_cyl_cfl_osh10_pe4.std.exo.2
#                                   slot_cyl_cfl_osh10_pe4.std.exo.3
#                      BIN_DIFF_PROG_ARGS -m
#                      BIN_RESULT out.e-s.0.4.0
#                                 out.e-s.0.4.1
#                                 out.e-s.0.4.2
#                                 out.e-s.0.4.3
#                      BIN_DIFF_PROG_CONF exodiff.cfg
#                      LABELS diagcg)
endif()

# Virtualization = 0.7

add_regression_test(cfl_u0.7 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_cfl.q unitsquare_01_3.6k.exo
                    ARGS -c slot_cyl_cfl.q -i unitsquare_01_3.6k.exo -v -u 0.7
                    BIN_BASELINE slot_cyl_cfl_pe1_u0.7.std.exo.0
                                 slot_cyl_cfl_pe1_u0.7.std.exo.1
                                 slot_cyl_cfl_pe1_u0.7.std.exo.2
                    BIN_RESULT out.e-s.0.3.0
                               out.e-s.0.3.1
                               out.e-s.0.3.2
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS diagcg)

# Parallel + Virtualization

add_regression_test(cfl_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl_cfl_diagcg.q unitsquare_01_3.6k.exo
                    ARGS -c slot_cyl_cfl_diagcg.q -i unitsquare_01_3.6k.exo
                         -v -u 0.5
                    BIN_BASELINE slot_cyl_cfl_diagcg_pe4_u0.5.std.exo.0
                                 slot_cyl_cfl_diagcg_pe4_u0.5.std.exo.1
                                 slot_cyl_cfl_diagcg_pe4_u0.5.std.exo.2
                                 slot_cyl_cfl_diagcg_pe4_u0.5.std.exo.3
                                 slot_cyl_cfl_diagcg_pe4_u0.5.std.exo.4
                                 slot_cyl_cfl_diagcg_pe4_u0.5.std.exo.5
                                 slot_cyl_cfl_diagcg_pe4_u0.5.std.exo.6
                                 slot_cyl_cfl_diagcg_pe4_u0.5.std.exo.7
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
                    LABELS diagcg)

# Parallel + Virtualization + Migration

add_regression_test(cfl_u0.9_migr ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl_cfl_diagcg.q unitsquare_01_3.6k.exo
                    ARGS -c slot_cyl_cfl_diagcg.q -i unitsquare_01_3.6k.exo
                         -v -l 10 -u 0.9 +balancer RandCentLB +LBDebug 1 +cs
                    BIN_BASELINE slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.0
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.1
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.2
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.3
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.4
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.5
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.6
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.7
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.8
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.9
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.10
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.11
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.12
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.13
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.14
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.15
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.16
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.17
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.18
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.19
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.20
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.21
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.22
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.23
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.24
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.25
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.26
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.27
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.28
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.29
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.30
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.31
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.32
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.33
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.34
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.35
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.36
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.37
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.38
                                 slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.39
                    BIN_RESULT out.e-s.0.40.0
                               out.e-s.0.40.1
                               out.e-s.0.40.2
                               out.e-s.0.40.3
                               out.e-s.0.40.4
                               out.e-s.0.40.5
                               out.e-s.0.40.6
                               out.e-s.0.40.7
                               out.e-s.0.40.8
                               out.e-s.0.40.9
                               out.e-s.0.40.10
                               out.e-s.0.40.11
                               out.e-s.0.40.12
                               out.e-s.0.40.13
                               out.e-s.0.40.14
                               out.e-s.0.40.15
                               out.e-s.0.40.16
                               out.e-s.0.40.17
                               out.e-s.0.40.18
                               out.e-s.0.40.19
                               out.e-s.0.40.20
                               out.e-s.0.40.21
                               out.e-s.0.40.22
                               out.e-s.0.40.23
                               out.e-s.0.40.24
                               out.e-s.0.40.25
                               out.e-s.0.40.26
                               out.e-s.0.40.27
                               out.e-s.0.40.28
                               out.e-s.0.40.29
                               out.e-s.0.40.30
                               out.e-s.0.40.31
                               out.e-s.0.40.32
                               out.e-s.0.40.33
                               out.e-s.0.40.34
                               out.e-s.0.40.35
                               out.e-s.0.40.36
                               out.e-s.0.40.37
                               out.e-s.0.40.38
                               out.e-s.0.40.39
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS diagcg migration)
