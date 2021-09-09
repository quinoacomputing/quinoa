# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Taylor-Green with ALE/Helmholtz"

inciter

  nstep 10
  ttyi 1
  cfl 0.5

  scheme alecg

  ale
    dvcfl 1.0
    mesh_velocity helmholtz
    vortmult 1.0  # remove all vorticity
    maxit 100
    tolerance 1.0e-3
    bc_dirichlet
      sideset 1 2 3 4 5 6 end
    end
    bc_sym
      sideset 1 2 3 4 5 6 end
    end
  end

  partitioning
   algorithm mj
  end

  compflow
    mesh filename "unitcube_1k.exo" end
    depvar u
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
