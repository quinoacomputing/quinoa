# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Initial uniform mesh refinement"

inciter

  nstep 1    # Max number of time steps
  cfl   0.2   # CFL coefficient
  ttyi 1      # TTY output interval

  scheme alecg

  partitioning
    algorithm mj
  end

  transport
    physics advection
    problem slot_cyl
#    bc_dirichlet
#      sideset 1 2 3 end
#    end
  end

  amr
    t0ref true
    dtref false
    dtfreq 5

    initial ic
    initial uniform_derefine
    initial ic
    initial uniform
    refvar c end
    error hessian

  end

  field_output
    interval 1
  end

end
