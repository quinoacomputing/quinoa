# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sod shock-tube"

inciter

  nsteps 25
  cfl 0.5
  ttyi 5  # TTY output interval
  scheme fv
  limiter vertexbasedp1

  partitioning
    algorithm mj
  end

  multimat

    physics euler
    problem sod_shocktube
    depvar u

    prelax 0

    flux hll

    nmat 2
    material
      eos stiffenedgas
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
    interval 1
    format scientific
    error l2
  end

  field_output
    interval 25
    var elem
      material_indicator "material_indicator_numerical"
      density "density_numerical"
      pressure "pressure_numerical"
      x-velocity "x-velocity_numerical"
    end
  end

end
