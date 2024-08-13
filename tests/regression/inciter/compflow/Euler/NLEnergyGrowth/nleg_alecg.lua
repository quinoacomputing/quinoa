inciter = {

title = "Euler equations computing nonlinear energy growth",

  term = 1.0,
  ttyi = 1,       -- TTY output interval
  cfl = 0.8,
  scheme = "alecg",

  partitioning = "mj",

  compflow = {
    physics = "euler",
    problem = "nl_energy_growth",
    alpha = 0.25,
    betax = 1.0,
    betay = 0.75,
    betaz = 0.5,
    r0 = 2.0,
    ce = -1.0,
    kappa = 0.8
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
    interval = 5,
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
