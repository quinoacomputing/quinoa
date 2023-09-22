# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Advection of 2D Gaussian hump"

inciter

  nstep 25   # Max number of time steps
  dt   1.0e-4 # Time step size
  ttyi 5     # TTY output interval
  scheme dgp2

  transport
    physics advection
    problem gauss_hump
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
    interval  5
    format    scientific
    error l2
    error linf
  end

  field_output
    interval 25
    var elem analytic C1 "c0_numerical" end
  end

end
