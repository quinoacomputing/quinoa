# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Advection of 2D Gaussian hump"

inciter

  nstep 10   # Max number of time steps
  dt   2.0e-3 # Time step size
  ttyi 1     # TTY output interval
  scheme dg

  transport
    physics advection
    problem gauss_hump
    ncomp 1

    bc_dirichlet
      sideset 1 2 3 4 5 6 end
    end
  end

  diagnostics
    interval  1
    format    scientific
    error l2
  end

  field_output
    interval 10
  end

end
