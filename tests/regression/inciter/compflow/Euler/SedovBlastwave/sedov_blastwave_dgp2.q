# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sedov blast wave"

inciter

  nstep 20   # Max number of time steps
  cfl 0.8
  ttyi 5      # TTY output interval
  scheme dgp2
  limiter vertexbasedp1

  compflow

    physics euler
    problem sedov_blastwave
    depvar u

    alpha 0.1
    beta 1.0
    p0 10.0

    material
      gamma 1.4 end
    end

    bc_sym
      sideset 1 2 end
    end
    bc_extrapolate
      sideset 3 end
    end

  end

  diagnostics
    interval  5
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
