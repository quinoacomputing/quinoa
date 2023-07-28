# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Advection of 2D Gaussian hump"

inciter

  nstep 50
  dt 2.0e-3
  ttyi 10

  scheme alecg

  partitioning
    algorithm mj
  end

  transport
    physics advection
    problem gauss_hump
    ncomp 1
    bc_sym
      sideset 1 end
    end
  end

  diagnostics
    interval  1
    format    scientific
    error l2
  end

  field_output
    interval 10
    var analytic C1 "c0_numerical" end
  end

end
