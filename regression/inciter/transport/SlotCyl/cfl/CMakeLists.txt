# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Virtualization = 0

add_regression_test(cfl ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_cfl.q unitsquare_01_3.6k.exo exodiff.cfg
                               slot_cyl_cfl_pe1_u0.0.std.exo
                    ARGS -c slot_cyl_cfl.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE slot_cyl_cfl_pe1_u0.0.std.exo
                    BIN_RESULT out.0
                    BIN_DIFF_PROG_CONF exodiff.cfg)

if (HAS_OMEGA_H)
  add_regression_test(cfl_osh ${INCITER_EXECUTABLE}
                      NUMPES 1
                      INPUTFILES slot_cyl_cfl.q box10.osh exodiff.cfg
                                 slot_cyl_cfl_osh10_pe1.std.exo
                      ARGS -c slot_cyl_cfl.q -i box10.osh -v
                      BIN_BASELINE slot_cyl_cfl_osh10_pe1.std.exo
                      BIN_RESULT out.0
                      BIN_DIFF_PROG_CONF exodiff.cfg)

  # parallel with nranks = nparts
  add_regression_test(cfl_osh ${INCITER_EXECUTABLE}
                      NUMPES 2
                      INPUTFILES slot_cyl_cfl.q box10_2.osh exodiff.cfg
                                 slot_cyl_cfl_osh10_pe2.std.exo.0
                                 slot_cyl_cfl_osh10_pe2.std.exo.1
                      ARGS -c slot_cyl_cfl.q -i box10_2.osh -v
                      BIN_BASELINE slot_cyl_cfl_osh10_pe2.std.exo.0
                                   slot_cyl_cfl_osh10_pe2.std.exo.1
                      BIN_RESULT out.0 out.1
                      BIN_DIFF_PROG_CONF exodiff.cfg)

#  # parallel with nranks > nparts
#  add_regression_test(cfl_osh ${INCITER_EXECUTABLE}
#                      NUMPES 4
#                      INPUTFILES slot_cyl_cfl.q box10_2.osh exodiff.cfg
#                                 slot_cyl_cfl_osh10_pe4.std.exo.0
#                                 slot_cyl_cfl_osh10_pe4.std.exo.1
#                                 slot_cyl_cfl_osh10_pe4.std.exo.2
#                                 slot_cyl_cfl_osh10_pe4.std.exo.3
#                      ARGS -c slot_cyl_cfl.q -i box10_2.osh -v
#                      BIN_BASELINE slot_cyl_cfl_osh10_pe4.std.exo.0
#                                   slot_cyl_cfl_osh10_pe4.std.exo.1
#                                   slot_cyl_cfl_osh10_pe4.std.exo.2
#                                   slot_cyl_cfl_osh10_pe4.std.exo.3
#                      BIN_DIFF_PROG_ARGS -m
#                      BIN_RESULT out.0 out.1 out.2 out.3
#                      BIN_DIFF_PROG_CONF exodiff.cfg)
endif()

add_regression_test(cfl ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl_cfl.q unitcube_01_31k.exo
                               exodiff.cfg
                               slot_cyl_cfl_pe4_u0.0.std.exo.0
                               slot_cyl_cfl_pe4_u0.0.std.exo.1
                               slot_cyl_cfl_pe4_u0.0.std.exo.2
                               slot_cyl_cfl_pe4_u0.0.std.exo.3
                    ARGS -c slot_cyl_cfl.q -i unitcube_01_31k.exo -v
                    BIN_BASELINE slot_cyl_cfl_pe4_u0.0.std.exo.0
                                 slot_cyl_cfl_pe4_u0.0.std.exo.1
                                 slot_cyl_cfl_pe4_u0.0.std.exo.2
                                 slot_cyl_cfl_pe4_u0.0.std.exo.3
                    BIN_RESULT out.0
                               out.1
                               out.2
                               out.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg)

# Virtualization = 0.7

add_regression_test(cfl_u0.7 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl_cfl.q unitsquare_01_3.6k.exo exodiff.cfg
                               slot_cyl_cfl_pe1_u0.7.std.exo.0
                               slot_cyl_cfl_pe1_u0.7.std.exo.1
                    ARGS -c slot_cyl_cfl.q -i unitsquare_01_3.6k.exo -v -u 0.7
                    BIN_BASELINE slot_cyl_cfl_pe1_u0.7.std.exo.0
                                 slot_cyl_cfl_pe1_u0.7.std.exo.1
                    BIN_RESULT out.0
                               out.1
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg)

# Parallel + Virtualization

add_regression_test(cfl_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl_cfl_diagcg.q unitsquare_01_3.6k.exo
                               exodiff.cfg
                               slot_cyl_cfl_diagcg_pe4_u0.5.std.exo.0
                               slot_cyl_cfl_diagcg_pe4_u0.5.std.exo.1
                               slot_cyl_cfl_diagcg_pe4_u0.5.std.exo.2
                               slot_cyl_cfl_diagcg_pe4_u0.5.std.exo.3
                               slot_cyl_cfl_diagcg_pe4_u0.5.std.exo.4
                               slot_cyl_cfl_diagcg_pe4_u0.5.std.exo.5
                               slot_cyl_cfl_diagcg_pe4_u0.5.std.exo.6
                               slot_cyl_cfl_diagcg_pe4_u0.5.std.exo.7
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

# Parallel + Virtualization + Migration

add_regression_test(cfl_u0.9_migr ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl_cfl_diagcg.q unitsquare_01_3.6k.exo
                               exodiff.cfg
                               slot_cyl_cfl_diagcg_pe4_u0.9.std.exo.0
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
                    ARGS -c slot_cyl_cfl_diagcg.q -i unitsquare_01_3.6k.exo
                         -v -u 0.9 +balancer RandCentLB +LBDebug 1 +cs
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
                    BIN_RESULT out.0 out.1 out.2 out.3 out.4 out.5 out.6 out.7
                               out.8 out.9 out.10 out.11 out.12 out.13 out.14
                               out.15 out.16 out.17 out.18 out.19 out.20 out.21
                               out.22 out.23 out.24 out.25 out.26 out.27 out.28
                               out.29 out.30 out.31 out.32 out.33 out.34 out.35
                               out.36 out.37 out.38 out.39
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    LABELS migration)


if (HAS_ROOT)

  # Virtualization = 0, RootMeshWriter check
  add_regression_test(cfl_root ${INCITER_EXECUTABLE}
                      NUMPES 1
                      INPUTFILES slot_cyl_cfl_root.q unitsquare_01_3.6k.exo
                                 exodiff.cfg  slot_cyl_cfl_pe1_u0.0.std.exo
	              ARGS -c slot_cyl_cfl_root.q -i unitsquare_01_3.6k.exo -v
		      FILECONV_INPUT out.0.root
		      FILECONV_RESULT out.0.exo
	              BIN_BASELINE slot_cyl_cfl_pe1_u0.0.std.exo
		      BIN_RESULT out.0.exo
		      BIN_DIFF_PROG_CONF exodiff.cfg)

  # Virtualization = 0, Multiprocessor check
  add_regression_test(cfl_root ${INCITER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES slot_cyl_cfl_root.q unitcube_01_31k.exo
	                         exodiff.cfg
	                         slot_cyl_cfl_pe4_u0.0.std.exo.0
	                         slot_cyl_cfl_pe4_u0.0.std.exo.1
	                         slot_cyl_cfl_pe4_u0.0.std.exo.2
		                 slot_cyl_cfl_pe4_u0.0.std.exo.3
		      ARGS -c slot_cyl_cfl_root.q -i unitcube_01_31k.exo -v
		      FILECONV_RESULT out.0.exo
				      out.1.exo
				      out.2.exo
				      out.3.exo
		      FILECONV_INPUT out.0.root
                                     out.1.root
                                     out.2.root
                                     out.3.root
		      BIN_BASELINE slot_cyl_cfl_pe4_u0.0.std.exo.0
			           slot_cyl_cfl_pe4_u0.0.std.exo.1
			           slot_cyl_cfl_pe4_u0.0.std.exo.2
			           slot_cyl_cfl_pe4_u0.0.std.exo.3
		      BIN_RESULT out.0.exo
			         out.1.exo
			         out.2.exo
			         out.3.exo
                      BIN_DIFF_PROG_ARGS -m
                      BIN_DIFF_PROG_CONF exodiff.cfg)

endif()
