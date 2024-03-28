# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Initial uniform mesh refinement"

inciter

  nstep 1     # Max number of time steps
  cfl   0.8   # CFL coefficient
  ttyi 1      # TTY output interval

  scheme alecg

  partitioning
    algorithm mj
  end

  transport
    physics advection
    problem slot_cyl
#    bc_dirichlet
#      sideset 1 2 3 end
#    end
  end

  amr
    t0ref true
    #initial uniform
    #initial ic
    initial coords
    initial coords

    coords
      x- 0.5
      #x+ 0.25 x- 0.75
    end
  end

  field_output
    interval 1
  end

end
