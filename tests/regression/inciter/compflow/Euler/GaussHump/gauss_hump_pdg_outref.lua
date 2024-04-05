inciter = {

  title = "Advection of 2D Gaussian hump",

  nstep = 10,  -- Max number of time steps
  cfl = 0.8,
  ttyi = 1,      -- TTY output interval
  scheme = "pdg",

  compflow = {
    physics = "euler",
    problem = "gauss_hump_compflow"
  },

  material = {
    {
      gamma = { 1.66666666666667 } -- =5/3 ratio of specific heats
    }
  },

  bc = {
    {
      symmetry = { 1 },
      dirichlet = { 2 },
      farfield = { 3 },
      pressure = 1.0,
      density = 1.0,
      velocity = { 0.0, 0.0, 0.0 }
    }
  },

  pref = {
    ndofmax = 10,
    tolref = 0.5
  },

  diagnostics = {
    interval = 5,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    refined = true,
    interval = 5,
    elemvar = {
      "analytic",
      "density",
      "x-velocity",
      "y-velocity",
      "z-velocity",
      "specific_total_energy",
      "pressure"
   },
    elemalias = {
      "",
      "density_numerical",
      "x-velocity_numerical",
      "y-velocity_numerical",
      "z-velocity_numerical",
      "specific_total_energy_numerical",
      "pressure_numerical"
   }
  }

}
