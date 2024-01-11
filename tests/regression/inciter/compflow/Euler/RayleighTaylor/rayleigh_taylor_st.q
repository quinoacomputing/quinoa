# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Euler equations computing the stationary Rayleigh-Taylor MMS problem"

inciter

  nstep 10
  #term 1.0
  ttyi 1       # TTY output interval
  cfl 0.5

  partitioning
    algorithm mj
  end

  compflow

    physics euler
    problem rayleigh_taylor
    alpha 1.0
    betax 1.0
    betay 1.0
    betaz 1.0
    p0 1.0
    r0 1.0
    kappa 0.0
    sysfct false

    material
      gamma 1.66666666666667 end # =5/3 ratio of specific heats
    end

    bc_dirichlet
      sideset 1 2 3 4 5 6 end
    end

  end

  field_output
    interval 1
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
    #error linf
  end

end
