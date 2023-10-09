# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Coupled Euler equations"

inciter

  nstep 1
  ttyi 1
  cfl 0.5

  partitioning
   algorithm rcb
  end

  scheme oversetfe

  compflow

    physics euler

    problem rayleigh_taylor
    alpha 1.0
    betax 1.0
    betay 1.0
    betaz 1.0
    p0 2.0
    r0 2.0
    kappa 1.0

    mesh
      # depvars are automatically assigned and can be referenced
      # downstream to request output variables
      filename "freestream_BGmesh_15k.exo" # depvar: 'a'
      filename "sphere_OSmesh_12k.exo" # depvar: 'b' ...
    end

    ic
      density  1.0 end
      velocity 0.0 0.0 0.0 end
      pressure 10.0 end
    end

    material
      gamma 1.66666666666667 end
    end

    bc_dirichlet
      mesh 1 end
      sideset 1 2 3 4 5 6 end
    end

    bc_farfield
      mesh 2 end
      sideset 101 end
      pressure 10.0
      density 1.0
      velocity 0.0 0.0 0.0 end
    end

    intergrid_boundary
      mesh 2 end
      sideset 102 end
    end

  end

  diagnostics
    interval  1
    format    scientific
    error l2
  end

  field_output
    var
      A1 A2 A3 A4 A5
      B1 B2 B3 B4 B5
    end
    interval 5
  end

end
