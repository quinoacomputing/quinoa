# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Advection of 2D Gaussian hump"

inciter

  nstep 200   # Max number of time steps
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
  end

  plotvar
    interval 50
  end

end
