# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Vortical flow"

inciter

  nstep 10   # Max number of time steps
  dt   1.0e-4 # Time step size
  ttyi 5      # TTY output interval

  scheme dg

  partitioning
    algorithm mj
  end

  reorder true

  compflow

    physics euler
    problem vortical_flow
    depvar u

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

  amr
    initial uniform
    refvar u end
    error jump
  end

  diagnostics
    interval  1
    format    scientific
    error l2
  end

  plotvar
    interval 2
  end

end
