# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Advection of cylinder"

inciter

  nstep 100   # Max number of time steps
  dt   1.0e-3 # Time step size
  ttyi 10     # TTY output interval
  ctau 1.0    # FCT mass diffusivity
  scheme dg

  transport
    physics advection
    problem cyl_advect
    ncomp 1
    depvar c

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
    interval  25
    format    scientific
    error l2
  end

  plotvar
    interval 50
  end

end
