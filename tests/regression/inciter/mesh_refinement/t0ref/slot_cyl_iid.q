# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Initial uniform mesh refinement"

inciter

  nstep 1    # Max number of time steps
  cfl   0.2   # CFL coefficient
  ttyi 1      # TTY output interval

  scheme diagcg

  fct true

  partitioning
    algorithm mj
  end

  transport
    physics advection
    problem slot_cyl
  end

  amr
    t0ref true
    dtref false
    dtfreq 5

    initial ic
    initial ic
    #initial ic
    #initial ic
    initial uniform_derefine
    refvar c end
    error hessian

  end

  field_output
    var
      node analytic C1 "c0_numerical"
    end
    interval 1
  end

end
