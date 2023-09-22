# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Initial uniform mesh refinement"

inciter

  nstep 1    # Max number of time steps
  cfl   0.2   # CFL coefficient
  ttyi 1      # TTY output interval

  scheme diagcg

  fct true

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
    dtref false
    dtfreq 5

    #initial uniform
    #initial uniform
    initial ic
    initial ic
    #initial ic
    #initial ic
    #initial coords
    #initial coords
    #initial coords
    refvar c end
    error jump
    tol_refine 0.8

#    coords
#      #x- 0.5
#      x+ 0.5 #x- 0.75
#      #y+ 0.25 y- 0.75
#    end

    #edgelist
    # 1 2 3
    #end

  end

  field_output
    interval 1
  end

end
