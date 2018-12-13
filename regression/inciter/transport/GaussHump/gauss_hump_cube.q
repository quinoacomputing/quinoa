# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Advection of 2D Gaussian hump"

inciter

  nstep 10   # Max number of time steps
  dt   2.0e-3 # Time step size
  ttyi 1     # TTY output interval
  ctau 1.0    # FCT mass diffusivity
  scheme dg

  transport
    physics advection
    problem gauss_hump
    ncomp 1
    depvar c

    bc_dirichlet
      sideset 1 2 3 4 5 6 end
    end
  end

  diagnostics
    interval  1
    format    scientific
    error l2
  end

  plotvar
    interval 10
  end

end
