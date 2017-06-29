# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Euler equations computing vortical flow"

inciter

  term 1.0
  ttyi 10       # TTY output interval
  cfl 0.8

  partitioning
   algorithm mj
  end

  compflow

    physics euler
    problem vortical_flow

    artvisc 0.05

    alpha 0.1
    beta 1.0
    p0 10.0

    material
      id 1
      gamma 1.66666666666667 # =5/3 ratio of specific heats
    end

    bc_dirichlet
      sideset 1 2 3 4 5 6 end
    end

  end

  plotvar
    interval 10
  end

  diagnostics
    interval  1
    format    scientific
    error l2
  end

end
