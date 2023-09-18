# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sod shock-tube"

inciter

  nstep 20
  cfl 0.3
  ttyi 10       # TTY output interval
  scheme pdg
  limiter vertexbasedp1

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
    ndofmax 10
    tolref 0
  end

  diagnostics
    interval  5
    format    scientific
    error l2
  end

  field_output
    interval 10
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
