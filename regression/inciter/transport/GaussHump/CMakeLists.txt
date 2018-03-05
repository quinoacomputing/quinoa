# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Virtualization = 0

add_regression_test(gauss_hump ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES gauss_hump.q unitsquare_01_3.6k.exo exodiff.cfg
                               gauss_hump.std.exo
                    ARGS -c gauss_hump.q -i unitsquare_01_3.6k.exo -v
                    BIN_BASELINE gauss_hump.std.exo
                    BIN_RESULT out.0
                    BIN_DIFF_PROG_CONF exodiff.cfg)
