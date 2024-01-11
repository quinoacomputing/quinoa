# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sod shock-tube, ALECG, ALE"

# This ALE config for this mesh yields close to Lagrangian mesh velocity and
# should run until t=0.2 without a problem.

inciter

  nstep 20    # Max number of time steps
  #term 0.2   # Max physical time
  ttyi 1      # TTY output interval

  cfl 0.9

  scheme alecg

  partitioning
    algorithm rcb
  end

  ale
    dvcfl 1.0
    mesh_velocity fluid
    smoother laplace
    vortmult 0.0      # no vorticity scaling in Laplace smoother
    maxit 20
    tolerance 2.0e-2
    mesh_motion 0 end # mesh moves only in x
    bc_dirichlet
      sideset 1 3 end
    end
  end

  compflow
    physics euler
    mesh filename "tube.exo" end
    problem sod_shocktube
    material
      gamma 1.4 end
    end
    bc_sym
      sideset 2 end
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
