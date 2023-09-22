# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Test GCL with ALECG"

inciter

  nstep 10    # Max number of time steps
  ttyi 1      # TTY output interval

  cfl 0.5

  partitioning
   algorithm mj
  end

  scheme alecg

  ale
    dvcfl 1.0
    mesh_velocity sine
  end

  compflow
    mesh filename "rectangle_01_1.5k.exo" end
    physics euler
    problem user_defined
    ic
      density  1.0 end
      velocity 0.0 0.0 0.0 end
      pressure 1.0 end
    end
    material
      gamma 1.4 end
    end
    bc_sym
      sideset 2 4 5 6 end
    end
  end

  field_output
    var
      density
      x-velocity
      y-velocity
      z-velocity
      specific_total_energy
      pressure
    end
    interval 10
  end

  diagnostics
    interval  1
    format    scientific
    error l2
  end

end
