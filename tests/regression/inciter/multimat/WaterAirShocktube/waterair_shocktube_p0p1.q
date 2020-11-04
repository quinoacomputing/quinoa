# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Water-air shock-tube"

inciter

  nstep 25
  dt 1.0e-8
  ttyi 10     # TTY output interval
  scheme p0p1
  limiter superbeep1

  partitioning
    algorithm mj
  end

  multimat

    physics veleq
    problem waterair_shocktube
    depvar u

    prelax 1
    prelax_timescale 0.25

    nmat 2
    material
      gamma 4.4 1.4 end # ratio of specific heats
      cv 951.36 717.5 end # specific heat at const volume
      pstiff 6.0e8 0.0 end # sg-eos stiffness parameter
    end

    bc_extrapolate
      sideset 1 3 end
    end
    bc_sym
      sideset 2 4 5 6 end
    end

  end

  diagnostics
    interval 1
    format    scientific
    error l2
  end

  field_output
    interval 25
  end

end
