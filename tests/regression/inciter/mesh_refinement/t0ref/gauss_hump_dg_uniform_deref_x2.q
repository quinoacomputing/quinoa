# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Advection of 2D Gaussian hump"

inciter

  nstep 10  # Max number of time steps
  dt   5.0e-4 # Time step size
  ttyi 1     # TTY output interval

  scheme dg

  partitioning
    algorithm mj
  end

  transport
    physics advection
    problem gauss_hump
    ncomp 1

    bc_extrapolate
      sideset 1 end
    end
    bc_inlet
      sideset 2 end
    end
    bc_outlet
      sideset 3 end
    end
  end

  amr
    t0ref true
    dtref false
    dtfreq 5
    initial uniform
    initial uniform_derefine
    initial uniform
    initial uniform_derefine
    initial uniform
    refvar c end
    error jump
  end

  diagnostics
    interval  2
    format    scientific
    error l2
  end

  field_output
    interval 2
    var elem analytic C1 "c0_numerical" end
  end

end
