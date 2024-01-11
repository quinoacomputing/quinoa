# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Taylor-Green with ALE, Laplacian smoothing, and vorticity scaling"

inciter

  nstep 10
  ttyi 1
  cfl 0.5

  scheme alecg

  ale
    dvcfl 1.0
    mesh_velocity fluid
    smoother laplace
    vortmult 1.0
    maxit 10
    tolerance 1.0
    bc_dirichlet
      sideset 1 2 3 4 5 6 end
    end
  end

  partitioning
   algorithm mj
  end

  compflow
    # if only a single filename is listen in a mesh ... end
    # block, it is not really a coupling, only a way to define the
    # input mesh filename in the control file, however if more
    # than one filename is given, the first one is considered the
    # background mesh (src) and the rest are destinations
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

  amr
    t0ref true
    initial uniform
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
