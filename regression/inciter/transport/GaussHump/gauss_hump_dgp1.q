# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Advection of 2D Gaussian hump"

inciter

  nstep 200   # Max number of time steps
  dt   2.0e-4 # Time step size
  ttyi 10     # TTY output interval
  ctau 1.0    # FCT mass diffusivity
  scheme dgp1

  transport
    physics advection
    problem gauss_hump
    ncomp 1
    depvar c

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
    interval  25
    format    scientific
    error l2
    error linf
  end

  plotvar
    interval 50
  end

end
