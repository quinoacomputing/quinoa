# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sod shock-tube"

inciter

  nstep 25   # Max number of time steps
  dt   1.0e-3 # Time step size
  ttyi 10     # TTY output interval
  scheme p0p1
  limiter superbeep1

  partitioning
    algorithm mj
  end

  multimat

    physics veleq
    problem sod_shocktube
    depvar u

    nmat 2
    material
      gamma 1.4 1.4 end # ratio of specific heats
    end

    bc_extrapolate
      sideset 1 3 end
    end
    bc_sym
      sideset 2 4 5 6 end
    end

  end

  diagnostics
    interval  1
    format    scientific
    error l2
  end

  plotvar
    interval 25
  end

end
