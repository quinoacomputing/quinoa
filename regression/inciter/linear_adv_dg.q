# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Zalesak's slotted cylinder"

inciter

  nstep 1     # Max number of time steps
  dt   0.001  # Time step size
  ttyi 1      # TTY output interval
  ctau 1.0    # FCT mass diffusivity

  discretization
    scheme dg
  end

  transport
    physics advdiff
    problem shear_diff
    ncomp 1
    depvar c
    diffusivity 3.0 2.0 1.0 end
    u0 10.0 end
    lambda 0.5 1.0 end
  end

  plotvar
    interval 1
  end

end
