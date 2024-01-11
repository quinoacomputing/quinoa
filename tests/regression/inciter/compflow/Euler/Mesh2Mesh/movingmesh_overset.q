# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Moving overset mesh"

inciter

  nstep 50
  ttyi 5
  cfl 0.75

  partitioning
   algorithm rcb
  end

  scheme oversetfe

  compflow

    physics euler
    problem user_defined

    mesh
      # depvars are automatically assigned and can be referenced
      # downstream to request output variables
      filename "freestream_BGmesh_15k.exo" velocity 0 0 0 end # depvar: 'a'
      filename "sphere_OSmesh_12k.exo" velocity 5 0 0 end # depvar: 'b' ...
    end

    ic
      density  1.0 end
      velocity 0.0 0.0 0.0 end
      pressure 1.0 end

      box
        xmin -0.5 xmax 0.525
        ymin -0.5 ymax 1
        zmin -0.5 zmax 1

        density 10.0
        velocity 0 0 0 end
        pressure 50.0
      end
    end

    material
      gamma 1.66666666666667 end
    end

    bc_dirichlet
      mesh 1 end
      sideset 4 6 end
    end

    bc_sym
      mesh 1 2 end
      sideset 1 2 3 5 102 end
    end

    bc_farfield
      mesh 2 end
      sideset 101 end
      pressure 1.0
      density 1.0
      velocity 0.0 0.0 0.0 end
    end

  end

  diagnostics
    interval  10
    format    scientific
    error l2
  end

  field_output
    var
      A1 A2 A3 A4 A5
      B1 B2 B3 B4 B5
    end
    interval 25
  end

  history_output
    interval 10
    point p1 0.5 0.25 0.25 end
  end

end
