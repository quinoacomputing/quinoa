# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Advection of 2D Gaussian hump"

inciter
  nstep 10  # Max number of time steps
  cfl 0.8
  ttyi 5      # TTY output interval
  scheme pdg

  compflow
    physics euler
    problem gauss_hump_compflow
    depvar u

    material
      gamma 1.66666666666667 end # =5/3 ratio of specific heats
    end

    bc_extrapolate
      sideset 1 end
    end
    bc_dirichlet
      sideset 2 end
    end
    bc_outlet
      farfield_pressure 1.0
      sideset 3 end
    end
  end

  pref
    ndofmax 10
    tolref 0.5
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
