inciter = {

  title = "Vortical flow",

  nstep = 10,   -- Max number of time steps
  dt = 1.0e-4, -- Time step size
  ttyi = 5,      -- TTY output interval
  scheme = "dg",

  partitioning = "mj",
  pelocal_reorder = true,

  compflow = {
    physics = "euler",
    problem = "vortical_flow",

    alpha = 0.1,
    beta = 1.0,
    p0 = 10.0
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

  amr = {
    t0ref = true,
    initial = {"uniform"},
    error = "jump"
  },

  diagnostics = {
    interval = 1,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 2,
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
