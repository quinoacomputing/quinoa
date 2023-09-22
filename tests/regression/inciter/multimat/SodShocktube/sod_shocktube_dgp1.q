# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sod shock-tube"

inciter

  nstep 25   # Max number of time steps
  cfl 0.8
  ttyi 10     # TTY output interval
  scheme dgp1
  limiter vertexbasedp1
  shock_detector_coeff 1.0
  limsol_projection false

  partitioning
    algorithm mj
  end

  multimat

    physics euler
    problem sod_shocktube

    prelax 0

    nmat 2
    material
      id 1 2 end
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

  field_output
    interval 25
    var elem
      F1 "volfrac1_numerical"
      F2 "volfrac2_numerical"
      density "density_numerical" # bulk density
      pressure "pressure_numerical" # bulk pressure
      specific_total_energy "total_energy_density_numerical" # bulk specific total energy
      x-velocity "x-velocity_numerical"
      y-velocity "y-velocity_numerical"
      z-velocity "z-velocity_numerical"
    end
  end

end
