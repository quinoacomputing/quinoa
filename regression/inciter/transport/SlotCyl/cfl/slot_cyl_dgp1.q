# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Zalesak's slotted cylinder"

inciter

  nstep 5     # Max number of time steps
  dt   5.0e-4    # Time step size
  ttyi 1      # TTY output interval
  ctau 1.0    # FCT mass diffusivity
  scheme dgp1

  transport
    depvar c
    physics advection
    problem slot_cyl
    bc_dirichlet
      sideset 1 end
    end
  end

  plotvar
    interval 5
  end

end
