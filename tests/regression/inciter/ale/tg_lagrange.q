# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Taylor-Green with pure Lagrangian mesh motion"

inciter

  nstep 10
  ttyi 1
  cfl 0.5

  scheme alecg

  ale
    dvcfl 1.0
    mesh_velocity fluid
  end

  partitioning
   algorithm mj
  end

  compflow
    mesh filename "unitcube_1k.exo" end
    physics euler
    problem taylor_green
    material
      gamma 1.4 end
    end
    bc_dirichlet
      sideset 1 2 3 4 5 6 end
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
    interval 1
    format   scientific
    error    l2
  end

end
