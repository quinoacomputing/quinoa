# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Interface advection"

inciter

  nstep 2    # Max number of time steps
  dt 2.5e-7   # Time step size
  ttyi 1      # TTY output interval
  scheme fv
  limiter vertexbasedp1

  multimat

    physics veleq
    problem interface_advection
    depvar u

    nmat 3
    material
      id 1 2 3 end
      gamma 1.4 1.4 1.4 end
      cv 83.33 717.5 717.5 end
    end

    bc_extrapolate
      sideset 1 2 3 4 5 6 end
    end

  end

  field_output
    interval 2
    var elem
      material_indicator
    end
  end

end
