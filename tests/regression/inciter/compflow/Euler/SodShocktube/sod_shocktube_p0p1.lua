inciter = {

  title = "Sod shock-tube",

  nstep = 100,   -- Max number of time steps
  dt = 2.0e-3,   -- Time step size
  ttyi = 10,     -- TTY output interval
  scheme = "p0p1",
  limiter = "superbeep1",

  compflow = {
    physics = "euler",
    problem = "sod_shocktube"
  },

  material = {
    {
      gamma = { 1.4 } -- ratio of specific heats
    }
  },

  bc = {
    {
      extrapolate = { 1, 3 },
      symmetry= { 2, 4, 5, 6 }
    }
  },

  diagnostics = {
    interval = 1,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 50,
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
