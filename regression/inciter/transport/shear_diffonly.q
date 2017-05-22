# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Dispersion from a point source in simple shear flow"

inciter

  nstep 10     # Max number of time steps
  t0   0.1     # Start time
  term 0.2     # Max time
  cfl 0.5
  ttyi 1      # TTY output interval

  transport
    physics advdiff
    problem shear_diff
    ncomp 1
    depvar c
    diffusivity 3.0 2.0 1.0 end
    u0 0.0 end
    lambda 0.0 0.0 end
  end

  plotvar
    interval 5
  end

end
