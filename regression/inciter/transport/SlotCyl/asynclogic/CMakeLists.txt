# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Virtualization = 0

add_regression_test(asynclogic ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v
                    TEXT_BASELINE slot_cyl.q
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic ${INCITER_EXECUTABLE}
                    NUMPES 3
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic ${INCITER_EXECUTABLE}
                    NUMPES 5
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic ${INCITER_EXECUTABLE}
                    NUMPES 6
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic ${INCITER_EXECUTABLE}
                    NUMPES 7
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic ${INCITER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

# Virtualization = 0.5

add_regression_test(asynclogic_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.5
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.5
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 3
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.5
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.5
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 5
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.5
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 6
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.5
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 7
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.5
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic_u0.5 ${INCITER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.5
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

# Virtualization = 0.9

add_regression_test(asynclogic_u0.9 ${INCITER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.9
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic_u0.9 ${INCITER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.9
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic_u0.9 ${INCITER_EXECUTABLE}
                    NUMPES 3
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.9
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic_u0.9 ${INCITER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.9
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic_u0.9 ${INCITER_EXECUTABLE}
                    NUMPES 5
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.9
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic_u0.9 ${INCITER_EXECUTABLE}
                    NUMPES 6
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.9
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic_u0.9 ${INCITER_EXECUTABLE}
                    NUMPES 7
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.9
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)

add_regression_test(asynclogic_u0.9 ${INCITER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES slot_cyl.q unitcube_01_31k.exo noop.ndiff.cfg
                    ARGS -c slot_cyl.q -i unitcube_01_31k.exo -v -u 0.9
                    TEXT_BASELINE unitcube_01_31k.exo
                    TEXT_RESULT diag
                    TEXT_DIFF_PROG_ARGS --trunc
                    TEXT_DIFF_PROG_CONF noop.ndiff.cfg)
