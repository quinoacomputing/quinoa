inciter = {

  title = "Advection of 2D Gaussian hump",

  nstep = 50,
  dt = 2.0e-3,
  ttyi = 10,
  scheme = "alecg",
  partitioning = "mj",

  transport = {
    problem = "gauss_hump",
    ncomp = 1
  },

  bc = {
    {
      symmetry = {1}
    }
  },

  physics = "advection",

  diagnostics = {
    interval = 1,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 10,
    nodevar = {"analytic", "C1"},
    nodealias = {"", "c0_numerical"}
  }

}
