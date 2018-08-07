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

add_regression_test(cfl_u0.5_migr ${INCITER_EXECUTABLE}
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
                         +balancer RandCentLB +LBDebug 2 +cs
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
