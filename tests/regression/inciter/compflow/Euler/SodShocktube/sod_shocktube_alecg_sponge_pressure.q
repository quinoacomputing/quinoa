# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sod shock-tube with a sponge-pressure BC"

inciter

  nstep 10    # Max number of time steps
  term 0.2    # Max physical time
  ttyi 1      # TTY output interval
  cfl 0.5

  scheme alecg

  partitioning
    algorithm mj
  end

  compflow
    physics euler
    problem sod_shocktube

    material
      gamma 1.4 end
    end

    bc_sym
      sideset 2 4 5 6 end
    end
    sponge
      sideset 2 end
      pressure 0.01 end    # reduce pressure gradient at sset 2 by 1%
    end
  end

  field_output
    interval 10000
    var
      density "density_numerical"
      x-velocity "x-velocity_numerical"
      y-velocity "y-velocity_numerical"
      z-velocity "z-velocity_numerical"
      specific_total_energy "specific_total_energy_numerical"
      pressure "pressure_numerical"
   end
  end

  diagnostics
    interval  1
    format    scientific
    error l2
  end

  history_output
    interval  1
    point p1 0.1 0.05 0.025 end
    point p2 0.9 0.05 0.025 end
  end

end
