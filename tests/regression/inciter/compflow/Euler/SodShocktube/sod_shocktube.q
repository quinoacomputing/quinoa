# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sod shock-tube"

inciter

  nstep 10    # Max number of time steps
  term 0.2    # Max physical time
  ttyi 1      # TTY output interval

  cfl 0.5
  #dt   1.0e-4 # Time step size

  scheme diagcg
  fcteps 1.0e-6

  partitioning
    algorithm mj
  end

  compflow
    physics euler
    problem sod_shocktube

    material
      gamma 1.4 end
    end

    bc_dirichlet
      sideset 1 3 end
    end

    bc_sym
      sideset 2 4 5 6 end
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

end
