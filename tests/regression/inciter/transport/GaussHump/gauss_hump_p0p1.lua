inciter = {

  title = "Advection of 2D Gaussian hump",

  nstep = 50,   -- Max number of time steps
  dt = 2.0e-3,  -- Time step size
  ttyi = 1,     -- TTY output interval
  scheme = "p0p1",
  limiter = "superbeep1",

  transport = {
    physics = "advection",
    problem = "gauss_hump",
    ncomp = 1
  },

  bc = {
    {
      extrapolate = {1},
      inlet = {2},
      outlet = {3}
    }
  },

  diagnostics = {
    interval = 2,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    elemvar = {"analytic", "C1"},
    elemalias = {"", "c0_numerical"},
    interval = 25
  }

}
