inciter = {

  title = "Advection of 2D Gaussian hump",

  nstep = 10,  -- Max number of time steps
  dt = 1.0e-3, -- Time step size
  ttyi = 1,     -- TTY output interval
  scheme = "dg",

  partitioning = "mj",
  pelocal_reorder = true,

  transport = {
    physics = "advection",
    problem = "gauss_hump",
    ncomp = 1
  }

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
   dtref = true,
   dtref_uniform = true,
   dtfreq = 5,
   error = "jump"
  },

  diagnostics = {
    interval = 2,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 1,
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
