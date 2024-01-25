# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Solid-fluid interface advection from Favrie and Gavrilyuk JCP 2009"

inciter

  nstep 10
  cfl 0.8
  ttyi 1  # TTY output interval
  scheme p0p1
  limiter vertexbasedp1

  partitioning
    algorithm mj
  end

  multimat

    physics euler
    problem user_defined

    flux laxfriedrichs

    prelax 0

    intsharp 0

    nmat 2
    # Copper
    material
      eos smallshearsolid
      id 1 end
      gamma 4.22 end # ratio of specific heats
      cv 3978.0 end # specific heat at const volume
      pstiff 342.0e8 end  # sg-eos stiffness parameter
      mu 9.2e10 end  # shear modulus
    end
    # Air
    material
      id 2 end
      gamma 1.4 end
      cv 717.5 end
    end

    ic
      # background (right-side conditions)
      materialid 2 end
      pressure 1.0e5 end
      temperature 300.0 end
      velocity 100.0 0.0 0.0 end

      # left-side conditions
      box
        materialid 1
        xmin -1e-10 xmax 0.4
        ymin -1.0 ymax 1.0
        zmin -1.0 zmax 1.0
        pressure 1.0e5
        temperature 300.0
        velocity 100.0 0.0 0.0 end
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
    interval 5
    var elem
      F1 "volfrac1"
      density "density" # bulk density
      pressure "pressure" # bulk pressure
      specific_total_energy "total_energy_density" # bulk specific total energy
      x-velocity "x-velocity"
      y-velocity "y-velocity"
      z-velocity "z-velocity"
    end
  end

end
