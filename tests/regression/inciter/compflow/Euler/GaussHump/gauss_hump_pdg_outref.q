# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Advection of 2D Gaussian hump"

inciter
  nstep 10  # Max number of time steps
  cfl 0.8
  ttyi 1      # TTY output interval
  scheme pdg

  compflow
    physics euler
    problem gauss_hump_compflow

    material
      gamma 1.66666666666667 end # =5/3 ratio of specific heats
    end

    bc_sym
      sideset 1 end
    end
    bc_dirichlet
      sideset 2 end
    end
    bc_farfield
      pressure 1.0
      density 1.0
      velocity 0.0 0.0 0.0 end
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

  field_output
    refined true
    interval 5
    var elem
      analytic
      density "density_numerical"
      x-velocity "x-velocity_numerical"
      y-velocity "y-velocity_numerical"
      z-velocity "z-velocity_numerical"
      specific_total_energy "specific_total_energy_numerical"
      pressure "pressure_numerical"
   end
  end

end
