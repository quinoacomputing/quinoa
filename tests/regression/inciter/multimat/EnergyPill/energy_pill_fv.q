# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Energy pill"

inciter

  nstep 50
  ttyi 5      # TTY output interval
  cfl 0.9

  scheme fv
  limiter vertexbasedp1

  partitioning
    algorithm mj
  end

  multimat
    physics energy_pill
    prelax 1
    prelax_timescale 0.25

    nmat 2

    material
      id 1 2 end
      gamma 1.4 1.4 end
      cv 717.5 717.5 end
    end

    # units: m,kg,s,K
    ic
      materialid 1 end
      temperature 300.0 end
      velocity 0.0 0.0 0.0 end
      pressure 77.85e+3 end

      meshblock
        blockid 2
        materialid 2
        mass 1024.0
        volume 0.64
        energy_content 9.0e+9
        initiate impulse #linear

        #linear
        #  point 0 0 0.75 end
        #  init_time 5e-5
        #  front_width 0.375
        #  velocity 8.2e+3
        #end
      end
    end

    bc_sym
      sideset 1 3 end
    end
  end

  diagnostics
    interval  5
    format    scientific
    error     l2
  end

  field_output
    interval 10
    sideset 1 end
    var elem
      density
      pressure
      x-velocity
      y-velocity
      z-velocity
      specific_total_energy
    end
  end

end
