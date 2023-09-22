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
    ncomp 2
    diffusivity 3.0 2.0 1.0 1.0 2.0 3.0 end
    u0 10.0 15.0 end
    lambda 0.5 1.0 0.75 0.25 end

    bc_dirichlet
      sideset 1 2 3 4 5 6 end
    end
  end

  diagnostics
    interval 3
    error l2
    error linf
  end

  field_output
    interval 5
    var analytic C1 "c0_numerical" C2 "c1_numerical" end
  end

end
