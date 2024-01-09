# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Zalesak's slotted cylinder"

inciter

  nstep 10    # Max number of time steps
  dt   0.001  # Time step size
  ttyi 1      # TTY output interval

  scheme oversetfe

  compflow
    physics euler
    problem user_defined
    material
      gamma 1.66666666666667 end
    end
    mesh filename "unitcube_1k.exo" end
    mesh filename "unitcube_1k.exo" end
    ic
      density  1.0 end
      velocity 0.0 0.0 0.0 end
      pressure 1.0 end
    end
  end

  field_output
    interval 5
  end

end
