# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Material interface advection"

inciter

  nstep 50 # Max number of time steps
  dt 2.5e-7
  ttyi 1    # TTY output interval
  scheme dg

  multimat

    physics euler
    problem interface_advection
    depvar u

    prelax 0

    nmat 3
    material
      id 1 2 3 end
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

  field_output
    interval 25
    var elem
	F1 "volfrac1_numerical"
	F2 "volfrac2_numerical"
	F3 "volfrac3_numerical"
	density "density_numerical" # bulk density
	x-velocity "x-velocity_numerical"
	y-velocity "y-velocity_numerical"
	z-velocity "z-velocity_numerical"
	pressure "pressure_numerical" # bulk presssure
	specific_total_energy "total_energy_density_numerical"  # bulk specific total energy
    end
  end

end
