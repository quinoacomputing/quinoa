# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Advection of 2D Gaussian hump"

inciter

  nstep 2000  # Max number of time steps
  dt   1.0e-3 # Time step size
  ttyi 50     # TTY output interval
  ctau 1.0    # FCT mass diffusivity

  discretization
    scheme dg
  end

  transport
    physics advection
    problem gauss_hump
    ncomp 1
    depvar c

    bc_sym
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

  plotvar
    interval 100
  end

end
