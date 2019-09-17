# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sod shock-tube"

inciter

  nstep 20
  cfl 0.3
  ttyi 10       # TTY output interval
  scheme pdg
  limiter superbeep1

  compflow

    physics euler
    problem sod_shocktube
    depvar u

    material
      gamma 1.4 end # ratio of specific heats
    end

    bc_extrapolate
      sideset 1 3 end
    end
    bc_sym
      sideset 2 4 5 6 end
    end

  end

  partitioning
    algorithm mj
  end

  pref
    ndofmax 4
    tolref 0
  end

  diagnostics
    interval  5
    format    scientific
    error l2
  end

  plotvar
    interval 5
  end

end
