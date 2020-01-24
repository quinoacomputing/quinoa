# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Serial, stationary

add_regression_test(compflow_euler_rayleigh_taylor_st ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES rayleigh_taylor_st.q unitcube_1k.exo
                    ARGS -c rayleigh_taylor_st.q -i unitcube_1k.exo -v
                    BIN_BASELINE rayleigh_taylor_st.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_st.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF rayleigh_taylor_diag.ndiff.cfg
                    LABELS diagcg)

# Parallel + no virtualization, stationary

add_regression_test(compflow_euler_rayleigh_taylor_st ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES rayleigh_taylor_st.q unitcube_1k.exo
                    ARGS -c rayleigh_taylor_st.q -i unitcube_1k.exo -v
                    BIN_BASELINE rayleigh_taylor_st_pe4.std.exo.0
                                 rayleigh_taylor_st_pe4.std.exo.1
                                 rayleigh_taylor_st_pe4.std.exo.2
                                 rayleigh_taylor_st_pe4.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag_st.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF rayleigh_taylor_diag.ndiff.cfg
                    LABELS diagcg)

# Serial, non-stationary

add_regression_test(compflow_euler_rayleigh_taylor ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES rayleigh_taylor.q unitcube_1k.exo
                    ARGS -c rayleigh_taylor.q -i unitcube_1k.exo -v
                    BIN_BASELINE rayleigh_taylor.std.exo
                    BIN_RESULT out.e-s.0.1.0
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF rayleigh_taylor_diag.ndiff.cfg
                    LABELS diagcg)

# Parallel + no virtualization, non-stationary

add_regression_test(compflow_euler_rayleigh_taylor ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES rayleigh_taylor.q unitcube_1k.exo
                    ARGS -c rayleigh_taylor.q -i unitcube_1k.exo -v
                    BIN_BASELINE rayleigh_taylor_pe4.std.exo.0
                                 rayleigh_taylor_pe4.std.exo.1
                                 rayleigh_taylor_pe4.std.exo.2
                                 rayleigh_taylor_pe4.std.exo.3
                    BIN_RESULT out.e-s.0.4.0
                               out.e-s.0.4.1
                               out.e-s.0.4.2
                               out.e-s.0.4.3
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_ARGS -m
                    BIN_DIFF_PROG_CONF exodiff.cfg
                    TEXT_BASELINE diag.std
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_CONF rayleigh_taylor_diag.ndiff.cfg
                    LABELS diagcg)
