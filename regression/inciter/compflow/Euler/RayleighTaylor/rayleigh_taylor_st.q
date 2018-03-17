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

    depvar c
    physics euler
    problem rayleigh_taylor
    alpha 1.0
    betax 1.0
    betay 1.0
    betaz 1.0
    p0 1.0
    r0 1.0
    kappa 0.0

    material
      id 1
      gamma 1.66666666666667 # =5/3 ratio of specific heats
    end

    bc_dirichlet
      sideset 1 2 3 4 5 6 end
    end

  end

  plotvar
    interval 1
  end

  diagnostics
    interval  1
    format    scientific
    error l2
    #error linf
  end

end
