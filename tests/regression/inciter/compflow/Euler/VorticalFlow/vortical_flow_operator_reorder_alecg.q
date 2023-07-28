# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Euler equations computing vortical flow"

inciter

  term 1.0
  ttyi 1       # TTY output interval
  cfl 0.8
  scheme alecg
  operator_reorder true

  partitioning
   algorithm mj
  end

  compflow

    physics euler
    problem vortical_flow

    alpha 0.1
    beta 1.0
    p0 10.0

    material
      gamma 1.66666666666667 end # =5/3 ratio of specific heats
    end

    bc_dirichlet
      sideset 1 2 3 4 5 6 end
    end

  end

  field_output
    interval 10
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
