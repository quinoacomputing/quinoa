# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Advection of cylinder"

inciter

  nstep 50    # Max number of time steps
  dt   1.0e-3 # Time step size
  ttyi 10     # TTY output interval
  scheme dgp1
  limiter superbeep1

  transport
    physics advection
    problem cyl_advect
    ncomp 1

    bc_extrapolate
      sideset 1 end
    end
    bc_dirichlet
      sideset 2 end
    end
    bc_outlet
      sideset 3 end
    end
  end

  diagnostics
    interval  10
    format    scientific
    error l2
  end

  field_output
    interval 50
    var elem analytic C1 "c0_numerical" end
  end

end
