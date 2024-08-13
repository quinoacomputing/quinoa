inciter = {

  title = "Sod shock-tube",

  nstep = 10,    -- Max number of time steps
  term = 0.2,    -- Max physical time
  ttyi = 1,      -- TTY output interval
  cfl = 0.5,
  scheme = "alecg",

  partitioning = "mj",

  compflow = {
    physics = "euler",
    problem = "sod_shocktube"
  },

  material = {
    {
      gamma = { 1.4 }
    }
  },

  bc = {
    {
      symmetry = { 2, 4, 5, 6 }
    }
  },

  field_output = {
    interval = 10000,
    time_interval = 1.0e-02,
    sideset = { 2, 4, 5, 6 },
    nodevar = {
      "density",
      "x-velocity",
      "y-velocity",
      "z-velocity",
      "specific_total_energy",
      "pressure"
   },
    nodealias = {
      "density_numerical",
      "x-velocity_numerical",
      "y-velocity_numerical",
      "z-velocity_numerical",
      "specific_total_energy_numerical",
      "pressure_numerical"
    }
  },

  diagnostics = {
    interval = 1,
    format = "scientific",
    error = "l2"
  },

}
