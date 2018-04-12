# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Zalesak's slotted cylinder"

inciter

  nstep 250   # Max number of time steps
  cfl   0.8   # CFL coefficient
  ttyi 10     # TTY output interval
  ctau 1.0    # FCT mass diffusivity

  discretization
    scheme diagcg
  end

  transport
    physics advection
    problem slot_cyl
  end

  plotvar
    interval 50
    filetype exodusii
  end

end
