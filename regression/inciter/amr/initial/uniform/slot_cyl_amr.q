# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Initial uniform mesh refinement"

inciter

  nstep 1     # Max number of time steps
  cfl   0.8   # CFL coefficient
  ttyi 1      # TTY output interval

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
