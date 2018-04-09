# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Vortical flow"

inciter

  nstep 100   # Max number of time steps
  dt   1.0e-4 # Time step size
  ttyi 5      # TTY output interval
  ctau 1.0    # FCT mass diffusivity

  discretization
    scheme dg
  end

  compflow

    physics euler
    problem vortical_flow

    alpha 0.1
    beta 1.0
    p0 10.0

    material
      id 1
      gamma 1.66666666666667 # =5/3 ratio of specific heats
    end

    bc_dirichlet
      sideset 1 2 3 end
    end

  end

  diagnostics
    interval  1
    format    scientific
    error l2
  end

  plotvar
    interval 10
  end

end
