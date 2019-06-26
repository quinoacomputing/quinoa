# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Material interface advection"

inciter

  nstep 50 # Max number of time steps
  dt 5e-7
  ttyi 1    # TTY output interval
  scheme dg

  multimat

    physics veleq
    problem interface_advection
    depvar u

    nmat 3
    material
      gamma 1.4 1.4 1.4 end
      cv 83.33 717.5 717.5 end
    end

    bc_sym
      sideset 1 end
    end
    bc_dirichlet
      sideset 2 end
    end
    bc_extrapolate
      sideset 3 end
    end

  end

  diagnostics
    interval  1 
    format    scientific
    error l2
  end

  plotvar
    interval 25
  end

end
