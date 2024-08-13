inciter = {

  title = "Zalesak's slotted cylinder",

  nstep = 10,    -- Max number of time steps
  dt = 0.001,    -- Time step size
  ttyi = 1,      -- TTY output interval
  scheme = "dg",

  transport = {
    problem = "slot_cyl"
  },

  physics = "advection",

  bc = {
    {
      dirichlet = {1, 2, 3, 4, 5, 6}
    }
  },

  field_output = {
    interval = 5
  }

}
