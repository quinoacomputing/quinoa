inciter = {

  title = "Euler equations computing vortical flow",

  term = 1.0,
  ttyi = 1,       -- TTY output interval
  cfl = 0.8,
  scheme = "alecg",
  operator_reorder = true,

  partitioning = "mj",

  compflow = {
    physics = "euler",
    problem = "vortical_flow",

    alpha = 0.1,
    beta = 1.0,
    p0 = 10.0,
  },

  material = {
    {
      gamma = { 1.66666666666667 } -- =5/3 ratio of specific heats
    }
  },

  bc = {
    {
      dirichlet = { 1, 2, 3, 4, 5, 6 }
    }
  },

  field_output = {
    interval = 10,
    nodevar = {
      "analytic",
      "density",
      "x-velocity",
      "y-velocity",
      "z-velocity",
      "specific_total_energy",
      "pressure"
    },
    nodealias = {
      "",
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
  }

}