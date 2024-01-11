# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Blast tube"

inciter

  nstep 100
  ttyi 5      # TTY output interval

  cfl 0.9

  scheme alecg

  partitioning
    algorithm mj
  end

  compflow
    physics euler

    # units: m,kg,s,
    ic
      density 0.912 end # kg/m^3
      velocity 0.0 0.0 0.0 end # m/s
      pressure 77.85e+3 end # Pa = N/m^2 = kg/m/s^2

      box
        xmin -0.200005 xmax 1.0e-14
        ymin -0.200005 ymax 0.200005
        zmin -0.000005 zmax 6.000005
        mass 768.48
        energy_content 9.0e+9

        initiate linear

        # linear propagation of energy source
        linear
          velocity -8.2e+3 # detonation velocity: 0.82 cm/us = 82e+3 m/s
          front_width 0.08
        end
      end
    end

    material
      gamma 1.4 end
      cv 717.5 end
    end

    bc_sym
      sideset 1 3 end
    end
  end

  field_output
    interval 25
    var
      density "density_numerical"
      specific_total_energy "specific_total_energy_numerical"
      pressure "pressure_numerical"
    end
  end

  diagnostics
    interval  5
    format    scientific
    error     l2
  end

end
