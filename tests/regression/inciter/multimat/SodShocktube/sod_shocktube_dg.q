# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Multi-material Sod shock-tube"

inciter

  nstep 100  # Max number of time steps
  dt   2.0e-4 # Time step size
  ttyi 10     # TTY output interval
  scheme dg

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
    interval 20
  end

end
