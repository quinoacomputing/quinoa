# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

foreach(scheme matcg diagcg dg alecg)
  foreach(virt 0.0 0.5 0.9)
    foreach(npes RANGE 1 8)
      add_regression_test(asynclogic_${scheme}_${virt} ${INCITER_EXECUTABLE}
        NUMPES ${npes}
        INPUTFILES slot_cyl_${scheme}.q unitcube_01_31k.exo noop.ndiff.cfg
        ARGS -c slot_cyl_${scheme}.q -i unitcube_01_31k.exo -v -u ${virt}
        TEXT_BASELINE noop.ndiff.cfg
        TEXT_RESULT diag
        TEXT_DIFF_PROG_ARGS --trunc
        TEXT_DIFF_PROG_CONF noop.ndiff.cfg)
    endforeach()
  endforeach()
endforeach()
