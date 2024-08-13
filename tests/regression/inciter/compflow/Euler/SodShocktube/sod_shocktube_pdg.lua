inciter = {

  title = "Sod shock-tube",

  nstep = 20,
  cfl = 0.3,
  ttyi = 10,       -- TTY output interval
  scheme = "pdg",
  limiter = "vertexbasedp1",

  partitioning = "mj",

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

  pref = {
    ndofmax = 10,
    tolref = 0
  },

  diagnostics = {
    interval = 5,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 10,
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
