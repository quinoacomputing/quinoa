# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Uniform mesh refinement during time stepping"

inciter

  nstep 9     # Max number of time steps
  cfl   0.8   # CFL coefficient
  ttyi 1      # TTY output interval

  scheme diagcg

  partitioning
    algorithm mj
  end

  transport
    physics advection
    problem slot_cyl
  end

  amr
    dtref true
    dtref_uniform true
    dtfreq 5
    refvar c end
    error jump
  end

  field_output
    interval 2
    var analytic C1 "c0_numerical" end
  end

end
