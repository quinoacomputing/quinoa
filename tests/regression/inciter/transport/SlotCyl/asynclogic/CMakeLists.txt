# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Add an asyncrhonous logic test for a discretization scheme. These tests are
# run in becnhmark mode, which means no large field output and thus they do not
# diff results. The main purpose of these tests is to see if the setup,
# initialization, and a few time steps go fine without a deadlock or some other
# problem. This can catch random or deterministic bugs in the asynchronous
# logic, especially if combined with Charm++'s randomized messsage queues as
# the tests are generated with a number of different combinations of PEs, degree
# of virtualization, and number of PEs per logical node (for Charm++'s SMP
# mode).
# Arguments:
#   scheme: which scheme to add a test for
#   virt: degree of virtualization: 0.0 ... 1.0 (both inclusive)
#   npes: total number of PEs to run test with
#   ppn: number of PEs per logical (compute) node, only used in Charm++'s
#        SMP mode, ignored in non-SMP mode
#   migr: 1: enable migration using RandCentLB, 0: no migration
function(add_async_test scheme virt npes ppn migr)
  if (migr)
    set(migr_cmd "+balancer RandCentLB +LBDebug 1 +cs")
    set(migr_append_name "_migr")
    set(migr_label "migration")
  endif()
  if(scheme STREQUAL "oversetfe")
    set(ARG "-c slot_cyl_${scheme}.q -v -u ${virt} -b ${migr_cmd}")
    set(TEXT_RES "diag.0;diag.1")
    set(TEXT_BAS "noop.ndiff.cfg;noop.ndiff.cfg")
  else()
    set(ARG "-c slot_cyl_${scheme}.q -i unitcube_1k.exo -v -u ${virt} -b ${migr_cmd}")
    set(TEXT_RES diag)
    set(TEXT_BAS "noop.ndiff.cfg")
  endif()
  add_regression_test(
    asynclogic_${scheme}_u${virt}${migr_append_name}
    ${INCITER_EXECUTABLE}
    NUMPES ${npes} PPN ${ppn}
    INPUTFILES slot_cyl_${scheme}.q unitcube_1k.exo
    ARGS ${ARG}
    TEXT_BASELINE ${TEXT_BAS}
    TEXT_RESULT ${TEXT_RES}
    TEXT_DIFF_PROG_ARGS --trunc
    TEXT_DIFF_PROG_CONF noop.ndiff.cfg
    LABELS ${scheme} ${migr_label})
endfunction()

foreach(scheme diagcg dg alecg fv oversetfe)
  foreach(virt 0.0 0.5 0.9)
    foreach(migr 0 1)

      # Add some tests across a range of PEs (both non-SMP & SMP)
      set(npes 1 3 8)
      foreach(n ${npes})
        add_async_test(${scheme} ${virt} ${n} 1 ${migr})
      endforeach()

      # In Charm++ SMP mode add some extra tests with npe != ppn
      if (CHARM_SMP)
        # Set npes & ppn as 4&2, 6&3, 6&2, 8&4, and 8&2, to correspond to
        # logical nodes x worker threads as 2x2, 2x3, 3x2, 2x4, and 4x2,
        # respectively. See also add_regresion_test.cmake for more details on
        # how to set NUMPES and PPN in SMP mode (npes must be a multiple of
        # ppn).
        set(npes 4 6 6 8 8)
        set(ppns 2 3 2 4 2)
        math(EXPR i "0")
        foreach(n ${npes})
          list(GET ppns ${i} p)
          add_async_test(${scheme} ${virt} ${n} ${p} ${migr})
          math(EXPR i "${i}+1")
        endforeach()
      endif()

    endforeach()
  endforeach()
endforeach()
