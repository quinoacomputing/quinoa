# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Time dependent BC test"

inciter

  term 0.5     # Max physical time, s
  ttyi 10      # TTY output interval
  cfl 0.9
  scheme oversetfe

  partitioning
    algorithm mj
  end

  compflow
    physics euler

    ic
      density 1.0 end
      velocity 0.0 0.0 0.0 end
      pressure 1.0 end
    end

    material
      gamma 1.4 end
      cv 717.5 end
    end

    bc_timedep
      sideset 1 end
      # The pressure, density, and velocity components are specified here as
      # a discrete function in time (tabular form). This table is sampled
      # to apply time dependent boundary conditions on the above side sets.
      fn # t    p    rho          u            v    w
        0.100   1.0  1.0          0.0          0.0  0.0
        0.101   5.0  2.818181818  1.606438658  0.0  0.0
        0.300   1.0  1.0          0.0          0.0  0.0
      end
    end

    bc_farfield
      pressure 1.0
      density 1.0
      velocity 0.0 0.0 0.0 end
      sideset 2 3 4 5 6 end
    end
  end

  field_output
    time_interval 0.1
    var node
      density
      pressure
      x-velocity
      y-velocity
      z-velocity
      specific_total_energy
    end
  end

  diagnostics
    interval  10
    format    scientific
    error     l2
  end

end
