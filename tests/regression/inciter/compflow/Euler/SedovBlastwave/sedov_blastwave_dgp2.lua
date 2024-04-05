inciter = {

  title = "Sedov blast wave",

  nstep = 20,   -- Max number of time steps
  cfl = 0.8,
  ttyi = 5,      -- TTY output interval
  scheme = "dgp1",
  limiter = "vertexbasedp1",

  compflow = {
    physics = "euler",
    problem = "sedov_blastwave"
  },

  material = {
    {
      gamma = { 1.4 }
    }
  },

  bc = {
    {
      symmetry = { 1, 2 },
      extrapolate = { 3 }
    }
  },

  diagnostics = {
    interval = 5,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 20,
    elemvar = {
      "density",
      "x-velocity",
      "y-velocity",
      "z-velocity",
      "specific_total_energy",
      "pressure"
   },
   elemalias = {
      "density_numerical",
      "x-velocity_numerical",
      "y-velocity_numerical",
      "z-velocity_numerical",
      "specific_total_energy_numerical",
      "pressure_numerical"
   }
  }

}
