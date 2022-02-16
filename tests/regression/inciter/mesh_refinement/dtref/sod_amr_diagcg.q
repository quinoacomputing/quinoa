# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sod AMR"

inciter

  nstep 1  # Max number of time steps
  dt 2e-5
  ttyi 1  # TTY output interval

  scheme diagcg
  fct true

  partitioning
    algorithm mj
  end

  compflow
    depvar u
    physics euler
    problem sod_shocktube

    material
      gamma 1.4 end
    end
  end

  amr
    maxlevels 1
    t0ref true
    dtref true
    dtfreq 1

    initial uniform
    refvar u end
    error jump
  end

  field_output
    interval 1
    var node
      density "density_numerical"
    end
  end

end
