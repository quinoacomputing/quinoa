inciter = {

  title = "Advection of 2D Gaussian hump",

  nstep = 10,  -- Max number of time steps
  dt = 5.0e-4, -- Time step size
  ttyi = 1,    -- TTY output interval
  scheme = "dg",

  partitioning = "mj",

  transport = {
    physics = "advection",
    problem = "gauss_hump",
    ncomp = 1
  },

  bc = {
    {
      extrapolate = { 1 },
      inlet = {
        sideset = { 2 }
      },
      outlet = { 3 }
    }
  },

  amr = {
    t0ref = true,
    dtref = false,
    dtfreq = 5,
    initial = {"uniform"},
    error = "jump"
  },

  diagnostics = {
    interval = 2,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 2,
    elemvar = {
      "analytic",
      "C1"
    },
    elemalias = {
      "",
      "c0_numerical"
    }
  }

}
