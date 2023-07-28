# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Advection of 2D Gaussian hump"

inciter

  nstep 50   # Max number of time steps
  dt   2.0e-3 # Time step size
  ttyi 1     # TTY output interval
  scheme p0p1
  limiter superbeep1

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

  diagnostics
    interval  2
    format    scientific
    error l2
  end

  field_output
    var elem analytic C1 "c0_numerical" end
    interval 25
  end

end
