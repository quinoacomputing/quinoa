# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Zalesak's slotted cylinder with initial uniform refinement"

inciter

  nstep 5     # Max number of time steps
  cfl   0.8   # CFL coefficient
  ttyi 1      # TTY output interval
  ctau 1.0    # FCT mass diffusivity

  amr
    initial uniform
  end

  transport
    physics advection
    problem slot_cyl
  end

  plotvar
    interval 1
  end

end
