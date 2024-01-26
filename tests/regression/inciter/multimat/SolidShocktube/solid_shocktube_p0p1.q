# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Solid-solid (Cu-Cu) shock tube"

inciter

  term 3e-5
  nstep 20
  cfl 0.8
  ttyi 5  # TTY output interval
  scheme p0p1
  limiter vertexbasedp1

  partitioning
    algorithm mj
  end

  multimat

    physics euler
    problem user_defined

    flux laxfriedrichs

    prelax 1

    intsharp 1

    nmat 2
    # Copper
    material
      eos smallshearsolid
      id 1 2 end
      gamma 4.22 4.22 end # ratio of specific heats
      cv 3978.0 3978.0 end # specific heat at const volume
      pstiff 342.0e8 342.0e8 end  # sg-eos stiffness parameter
      mu 9.2e10 9.2e10 end  # shear modulus
    end

    ic
      # background (right-side conditions)
      materialid 2 end
      pressure 1e5 end
      temperature 343.85505 end
      velocity 0.0 0.0 0.0 end

      # left-side conditions
      box
        materialid 1
        xmin -1e-10 xmax 0.5
        ymin -1.0 ymax 1.0
        zmin -1.0 zmax 1.0
        pressure 5e9
        temperature 343.85505
        velocity 0.0 0.0 0.0 end
      end
    end

    bc_extrapolate
      sideset 1 3 end
    end
    bc_sym
      sideset 2 4 5 6 end
    end

  end

  diagnostics
    interval 5
    format    scientific
    error l2
  end

  field_output
    interval 10
    var elem
      F1
      density # bulk density
      pressure # bulk pressure
      specific_total_energy # bulk specific total energy
      x-velocity
      y-velocity
      z-velocity
    end
  end

end
