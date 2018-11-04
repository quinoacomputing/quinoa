# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sod shock-tube"

inciter

  nstep 100   # Max number of time steps
  dt   2.0e-3 # Time step size
  ttyi 10     # TTY output interval
  scheme dg

  compflow

    physics euler
    problem sod_shocktube
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
    interval  1
    format    scientific
    error l2
  end

  plotvar
    interval 20
  end

end
