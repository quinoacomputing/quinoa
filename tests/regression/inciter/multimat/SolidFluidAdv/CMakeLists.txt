# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Parallel
# add_regression_test(multimat_solidfluid_adv_p0p1 ${INCITER_EXECUTABLE}
#                     NUMPES 4
#                     INPUTFILES solidfluid_adv_p0p1.lua rectangle_01_1.5k.exo
#                     ARGS -c solidfluid_adv_p0p1.lua -i rectangle_01_1.5k.exo -v
#                     BIN_BASELINE solidfluid_adv_p0p1_pe4.std.exo.0
#                                  solidfluid_adv_p0p1_pe4.std.exo.1
#                                  solidfluid_adv_p0p1_pe4.std.exo.2
#                                  solidfluid_adv_p0p1_pe4.std.exo.3
#                     BIN_RESULT out.e-s.0.4.0
#                                out.e-s.0.4.1
#                                out.e-s.0.4.2
#                                out.e-s.0.4.3
#                     BIN_DIFF_PROG_CONF exodiff_dg.cfg
#                     BIN_DIFF_PROG_ARGS -m
#                     TEXT_BASELINE diag_p0p1.std
#                     TEXT_RESULT diag
#                     TEXT_DIFF_PROG_CONF ndiff.cfg
#                     LABELS dg multimat solids)
