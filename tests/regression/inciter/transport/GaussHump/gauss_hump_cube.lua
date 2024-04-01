inciter = {

  title = "Advection of 2D Gaussian hump",

  nstep = 10,   -- Max number of time steps
  dt = 2.0e-3,  -- Time step size
  ttyi = 1,     -- TTY output interval
  scheme = "dg",

  transport = {
    problem = "gauss_hump",
    ncomp = 1
  },

  physics = "advection",

  bc = {
    {
      dirichlet = {1, 2, 3, 4, 5, 6}
    }
  },

  diagnostics = {
    interval = 1,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 10
  }

}
