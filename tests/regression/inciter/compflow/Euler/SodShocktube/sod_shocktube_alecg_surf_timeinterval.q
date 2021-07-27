# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sod shock-tube"

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
    depvar u
    physics euler
    problem sod_shocktube

    material
      gamma 1.4 end
    end

    bc_sym
      sideset 2 4 5 6 end
    end
  end

  field_output
    interval 10000
    time_interval 1.0e-02
    sideset 2 4 5 6 end
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

end
