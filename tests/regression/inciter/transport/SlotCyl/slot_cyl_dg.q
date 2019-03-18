# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Zalesak's slotted cylinder"

inciter

  nstep 5       # Max number of time steps
  dt   5.0e-4   # Time step size
  ttyi 1        # TTY output interval
  scheme dg

  transport
    depvar c
    physics advection
    problem slot_cyl
    bc_dirichlet
      sideset 1 end
    end
  end

  plotvar
    interval 1
  end

end
