# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Zalesak's slotted cylinder"

inciter

  nstep 10    # Max number of time steps
  dt   0.001  # Time step size
  ttyi 1      # TTY output interval
  ctau 1.0    # FCT mass diffusivity

  scheme alecg

  transport
    physics advection
    problem slot_cyl
  end

  field_output
    interval 5
  end

end
