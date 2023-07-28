# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sod shock-tube, ALECG, sine mesh motion"

inciter

  nstep 10    # Max number of time steps
  ttyi 1      # TTY output interval

  cfl 0.5

  scheme alecg

  partitioning
    algorithm mj
  end

  ale
    dvcfl 1.0
    mesh_velocity sine
  end

  compflow
    physics euler
    mesh filename "rectangle_01_1.5k.exo" end
    problem sod_shocktube
    material
      gamma 1.4 end
    end
    bc_sym
      sideset 2 4 5 6 end
    end
  end

  field_output
    interval 10
    var node
      density
      x-velocity
      y-velocity
      z-velocity
      specific_total_energy
      pressure
    end
  end

  diagnostics
    interval  1
    format    scientific
    error l2
  end

end
