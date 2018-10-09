# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Initial uniform mesh refinement"

inciter

  nstep 1     # Max number of time steps
  cfl   0.8   # CFL coefficient
  ttyi 1      # TTY output interval

  scheme diagcg

  partitioning
    algorithm mj
  end

  transport
    depvar c
    physics advection
    problem slot_cyl
#    bc_dirichlet
#      sideset 1 2 3 end
#    end
  end

  amr
    #initial uniform
    #initial uniform
    #initial ic
    #initial uniform
    #initial ic
    #initial uniform
    #initial ic
    initial coords
    initial coords
    #initial coords
    #initial ic
    #initial uniform
    #initial ic
    #initial uniform
    #initial ic
    #initial uniform
    #refvar c end
    #error jump

    coordref
      x- 0.5
      #x+ 0.25 x- 0.75
    end
    #initref
    # 1 2 3
    #end
  end

  plotvar
    interval 1
  end

end
