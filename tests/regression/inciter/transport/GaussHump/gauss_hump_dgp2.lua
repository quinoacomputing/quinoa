inciter = {

  title = "Advection of 2D Gaussian hump",

  nstep = 25,   -- Max number of time steps
  dt = 1.0e-4,  -- Time step size
  ttyi = 5,     -- TTY output interval
  scheme = "dgp2",

  transport = {
    problem = "gauss_hump",
    ncomp = 1
  },

  physics = "advection",

  bc = {
    {
      extrapolate = {1},
      dirichlet = {2},
      outlet = {3}
    }
  },

  diagnostics = {
    interval = 5,
    format = "scientific",
    error = "linf"
  },

  field_output = {
    interval = 25,
    elemvar = {"analytic", "C1"},
    elemalias = {"", "c0_numerical"}
  }

}
