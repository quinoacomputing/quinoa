# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Zalesak's slotted cylinder"

inciter

  nstep 50    # Max number of time steps
  cfl   0.8   # CFL coefficient
  ttyi 1      # TTY output interval
  ctau 1.0    # FCT mass diffusivity

  scheme diagcg

  transport
    physics advection
    problem slot_cyl
  end

  field_output
    interval 50
    filetype exodusii
    var analytic C1 "c0_numerical" end
  end

end
