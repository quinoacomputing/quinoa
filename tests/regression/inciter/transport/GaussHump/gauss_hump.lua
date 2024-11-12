inciter = {

  title = "Advection of 2D Gaussian hump",

  nstep = 100,  -- Max number of time steps
  dt = 2.0e-3,  -- Time step size
  ttyi = 1,     -- TTY output interval
  scheme = "dg",

  transport= {
    physics = "advection",
    problem = "gauss_hump",
    ncomp = 1
  },

  bc = {
    {
      extrapolate = {1},
      inlet = {
        sideset = { 2 }
      },
      outlet = {3}
    }
  },

  diagnostics = {
    interval = 2,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 10,
    elemvar = {"analytic", "C1"},
    elemalias = {"", "c0_numerical"}
  }

}
