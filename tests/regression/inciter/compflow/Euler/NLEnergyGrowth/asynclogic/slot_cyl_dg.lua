inciter = {

  title = "Euler equations computing nonlinear energy growth",

  nstep = 10,    -- Max number of time steps
  dt = 0.001,    -- Time step size
  ttyi = 1,      -- TTY output interval
  scheme = "dg",

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
      dirichlet = {1, 2, 3, 4, 5, 6}
    }
  },

  field_output = {
    interval = 5
  }

}
