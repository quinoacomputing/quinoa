# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Euler equations computing nonlinear energy growth"

inciter

  term 1.0
  ttyi 1       # TTY output interval
  cfl 0.8

  partitioning
    algorithm mj
  end

  scheme alecg

  compflow

    physics euler
    problem nl_energy_growth
    alpha 0.25
    betax 1.0
    betay 0.75
    betaz 0.5
    r0 2.0
    ce -1.0
    kappa 0.8

    material
      gamma 1.66666666666667 end # =5/3 ratio of specific heats
    end

    bc_dirichlet
      sideset 1 2 3 4 5 6 end
    end

  end

  field_output
    interval 5
    var
      analytic
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
