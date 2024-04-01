inciter = {

  title = "Interface advection",

  nstep = 2,     -- Max number of time steps
  dt = 2.5e-7,   -- Time step size
  ttyi = 1,      -- TTY output interval
  scheme = "fv",
  limiter = "vertexbasedp1",

  multimat = {
    problem = "interface_advection",
    nmat = 3
  },

  physics = "euler",

  material = {
    {
      id = { 1, 2, 3 },
      gamma = { 1.4, 1.4, 1.4 },
      cv = { 83.33, 717.5, 717.5 }
    }
  },

  bc = {
    {
      extrapolate = {1, 2, 3, 4, 5, 6}
    }
  },

  field_output = {
    interval = 2,
    elemvar = { "material_indicator" }
  }

}
