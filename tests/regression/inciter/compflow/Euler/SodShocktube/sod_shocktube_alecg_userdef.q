# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sod shock-tube"

inciter

  nstep 10    # Max number of time steps
  term 0.2    # Max physical time
  ttyi 1      # TTY output interval
  cfl 0.5

  scheme alecg

  partitioning
    algorithm mj
  end

  compflow
    physics euler

    #problem sod_shocktube
    ic
      density -1.0 end                  # overwritten by boxes
      velocity 100.0 100.0 100.0 end    # overwritten by boxes
      pressure -1.0 end                 # overwritten by boxes
      box
        xmin -0.5 xmax 0.5
        ymin -0.5 ymax 0.5
        zmin -0.5 zmax 0.5
        density 1.0
        pressure 1.0
      end
      box
        xmin  0.5 xmax 1.5
        ymin -0.5 ymax 0.5
        zmin -0.5 zmax 0.5
        density 0.125
        pressure 0.1
      end
    end

    material
      gamma 1.4 end
    end

    bc_sym
      sideset 2 4 5 6 end
    end
  end

  field_output
    interval 10000
    var
      density "density_numerical"
      x-velocity "x-velocity_numerical"
      y-velocity "y-velocity_numerical"
      z-velocity "z-velocity_numerical"
      specific_total_energy "specific_total_energy_numerical"
      pressure "pressure_numerical"
   end
  end

  diagnostics
    interval  1
    format    scientific
    error l2
  end

end
