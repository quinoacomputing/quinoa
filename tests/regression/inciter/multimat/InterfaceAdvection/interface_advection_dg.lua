inciter = {

  title = "Material interface advection",

  nstep = 50,  -- Max number of time steps
  dt = 2.5e-7,
  ttyi = 1,    -- TTY output interval
  scheme = "dgp0",
  lowspeed_kp = 1.0,

  multimat = {
    physics = "euler",
    problem = "interface_advection",
    prelax = 0,
    nmat = 3
  },

  material = {
    {
      id = { 1, 2, 3 },
      gamma = { 1.4, 1.4, 1.4 },
      cv = { 83.33, 717.5, 717.5 }
    }
  },

  bc = {
    {
      symmetry = { 1 },
      dirichlet = { 2 },
      extrapolate = { 3 }
    }
  },

  diagnostics = {
    interval = 1,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 25,
    elemvar = {
	    "F1",
	    "F2",
	    "F3",
	    "density",
	    "x-velocity",
	    "y-velocity",
	    "z-velocity",
	    "pressure",
	    "specific_total_energy"
    },
    elemalias = {
	    "volfrac1_numerical",
	    "volfrac2_numerical",
	    "volfrac3_numerical",
	    "density_numerical",
	    "x-velocity_numerical",
	    "y-velocity_numerical",
	    "z-velocity_numerical",
	    "pressure_numerical",
	    "total_energy_density_numerical"
    }
  }

}
