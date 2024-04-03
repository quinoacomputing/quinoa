# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set(multimat_solid_shocktube_p0p1_BIN_DIFF_PROG_CONF exodiff_dg_intel.cfg)
  set(multimat_solid_shocktube_p0p1_TXT_DIFF_PROG_CONF ndiff_intel.cfg)
else()
  set(multimat_solid_shocktube_p0p1_BIN_DIFF_PROG_CONF exodiff_dg.cfg)
  set(multimat_solid_shocktube_p0p1_TXT_DIFF_PROG_CONF ndiff.cfg)
endif()
add_regression_test(multimat_solid_shocktube_p0p1 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES solid_shocktube_p0p1.lua rectangle_01_int50_struc896.exo
                    ARGS -c solid_shocktube_p0p1.lua -i rectangle_01_int50_struc896.exo -v
                    BIN_BASELINE solid_shocktube_p0p1.std.exo.0
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF
                      ${multimat_solid_shocktube_p0p1_BIN_DIFF_PROG_CONF}
                    BIN_DIFF_PROG_ARGS -m
                    TEXT_BASELINE diag_p0p1.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF
                      ${multimat_solid_shocktube_p0p1_TXT_DIFF_PROG_CONF}
                    LABELS dg multimat solids)
