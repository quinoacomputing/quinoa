inciter = {

  title = "Vortical flow",

  nstep = 50,    -- Max number of time steps
  dt = 1.0e-5,   -- Time step size
  ttyi = 5,      -- TTY output interval
  scheme = "dgp1",
  flux = "laxfriedrichs",

  compflow = {
    physics = "euler",
    problem = "vortical_flow",
    alpha = 0.1,
    beta = 1.0,
    p0 = 10.0
  },

  material = {
    {
      gamma = { 1.66666666666667 }  -- =5/3 ratio of specific heats
    }
  },

  bc = {
    {
      dirichlet = { 1, 2, 3, 4, 5, 6 }
    }
  },

  diagnostics = {
    interval = 1,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    elemvar = {
      "analytic",
      "density",
      "x-velocity",
      "y-velocity",
      "z-velocity",
      "specific_total_energy",
      "pressure",
    },
    elemalias = {
      "",
      "density_numerical",
      "x-velocity_numerical",
      "y-velocity_numerical",
      "z-velocity_numerical",
      "specific_total_energy_numerical",
      "pressure_numerical"
    },
    interval = 10
  }

}
