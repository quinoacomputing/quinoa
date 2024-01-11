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
    problem rotated_sod_shocktube

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

  diagnostics
    interval  1
    format    scientific
    error l2
  end

  field_output
    interval 20
    var elem
      density "density_numerical"
      x-velocity "x-velocity_numerical"
      y-velocity "y-velocity_numerical"
      z-velocity "z-velocity_numerical"
      specific_total_energy "specific_total_energy_numerical"
      pressure "pressure_numerical"
    end
  end

end
