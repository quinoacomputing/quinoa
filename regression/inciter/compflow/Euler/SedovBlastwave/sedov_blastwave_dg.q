# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sedov blast-wave"

inciter

  nstep 400   # Max number of time steps
  dt   5.0e-4 # Time step size
  ttyi 10     # TTY output interval
  scheme dg

  compflow

    physics euler
    problem sedov_blastwave
    depvar u

    material
      id 1
      gamma 1.4 # ratio of specific heats
    end

    bc_extrapolate
      sideset 1 3 end
    end
    bc_sym
      sideset 2 4 5 6 end
    end

  end

  diagnostics
    interval  50
    format    scientific
    error l2
  end

  plotvar
    interval 50
  end

end
