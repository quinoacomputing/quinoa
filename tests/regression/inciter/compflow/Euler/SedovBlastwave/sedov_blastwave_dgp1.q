# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sedov blast wave"

inciter

  nstep 20   # Max number of time steps
  cfl 0.3
  ttyi 5      # TTY output interval
  scheme dgp1
  limiter superbeep1

  compflow

    physics euler
    problem sedov_blastwave
    depvar u

    alpha 0.1
    beta 1.0
    p0 10.0

    material
      id 1
      gamma 1.4
    end

    bc_sym
      sideset 1 2 end
    end
    bc_extrapolate
      sideset 3 end
    end

  end

  diagnostics
    interval  5
    format    scientific
    error l2
  end

  plotvar
    interval 20
  end

end
