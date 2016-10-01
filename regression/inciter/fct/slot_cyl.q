# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Zalesak's slotted cylinder"

inciter

  nstep 5     # Max number of time steps
  dt   0.001  # Time step size
  ttyi 1      # TTY output interval
  ctau 1.0    # FCT mass diffusivity

  transport
    physics advection
    problem slot_cyl
  end

  plotvar
    interval 1
  end

end
