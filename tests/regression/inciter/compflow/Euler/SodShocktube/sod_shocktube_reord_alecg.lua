inciter = {

  title = "Sod shock-tube",

  nstep = 10,    -- Max number of time steps
  term = 0.2,    -- Max physical time
  ttyi = 1,      -- TTY output interval
  cfl = 0.5,
  scheme = "alecg",
  operator_reorder = true,

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

  history_output = {
    interval = 1,
    point = {
      { id = "p1", coord = {0.1, 0.05, 0.025} },
      { id = "p2", coord = {0.9, 0.05, 0.025} }
    }
  }

}
